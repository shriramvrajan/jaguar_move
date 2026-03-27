#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector propagate_cpp(NumericVector attract_vec,  // flattened attract matrix (col-major), length slice_size * ncol_dest
                            IntegerVector to_dest_vec,  // flattened connectivity map (col-major, 1-indexed from R), same length
                            int ncell,             // cells per local neighborhood: (2*max_dist+1)^2
                            int nobs,              // number of observations in individual trajectory
                            int nsteps,            // propagation steps
                            double bg_rate,        // background movement rate for long steps
                            int ncol_dest) {       // columns in dest/to_dest matrix (= ncol of inner nbhd)

  // slice_size is the number of elements in one time-slice of the 3D array.
  // Each slice is a matrix of ncell rows (cells in local neighborhood) by nobs columns (observations).
  int slice_size = ncell * nobs;

  // Full output, flat vector representing 3D array current[cell, obs, step]
  //   current[c, o, s] = current[c + o*ncell + s*slice_size]
  int total_size = slice_size * nsteps;
  NumericVector current(total_size, 0.0);

  // Initialize step 0: probability 1.0 at center cell of each observation's neighborhood.
  // ncell always odd ((2k+1)^2), so integer division gives 0-based center index.
  int center = ncell / 2;
  for (int obs = 0; obs < nobs; obs++) {
    current[center + obs * ncell] = 1.0;
  }

  // Pre-allocate working memory. Reused every iteration — no allocation inside the loop.
  // p_step[i] will hold the total incoming probability to cell i (across all neighbor directions).
  std::vector<double> p_step(slice_size);
  // One normalization total per observation column.
  std::vector<double> col_totals(nobs);

  // Main propagation loop: steps 0 to nsteps-2 produce slices 1 to nsteps-1
  for (int j = 0; j < nsteps - 1; j++) {

    // Byte offsets into `current` for the source (j) and destination (j+1) slices.
    int offset_in  = j * slice_size;
    int offset_out = (j + 1) * slice_size;

    // Zero out the accumulator for this iteration.
    std::fill(p_step.begin(), p_step.end(), 0.0);

    // --- Fused gather + rowSum ---
    // For each destination cell (row of to_dest),
    // walk across its ncol_dest incoming connections, compute the source
    // probability * transition weight, and accumulate probability mass.
    for (int row = 0; row < slice_size; row++) {
      double incoming = 0.0;

      for (int c = 0; c < ncol_dest; c++) {
        // Value of to_dest[row, c] (1-indexed from R). If NA, skip this connection.
        int td_val = to_dest_vec[row + c * slice_size];
        if (td_val == NA_INTEGER) continue;

        // Convert from R 1-indexing to C++ 0-indexing.
        int src_idx = td_val - 1;

        // attract_vec may contain NA (represented as NaN) for out-of-bounds cells.
        // In R, these get dropped by na.rm = TRUE in rowSums. We skip them here.
        double attract_val = attract_vec[src_idx];
        if (std::isnan(attract_val)) continue;

        // The key computation. In R, step_prob was formed by recycling current[,,j]
        // (length slice_size) across attract[] (length slice_size * ncol_dest).
        // step_prob[k] = current[k % slice_size, , j] * attract[k]
        // So the source cell's current probability is at index (src_idx % slice_size)
        // within the j-th slice.
        int current_row = src_idx % slice_size;
        incoming += current[offset_in + current_row] * attract_val;
      }

      p_step[row] = incoming;
    }

    // Background rate and normalization
    // p_bg = p + bg - p*bg
    // This ensures every cell retains at least bg_rate probability of being occupied,
    // even if no directed movement reaches it. When bg_rate ~= 0, nearly no-op.

    // Apply bg_rate and compute column totals.
    std::fill(col_totals.begin(), col_totals.end(), 0.0);
    for (int obs = 0; obs < nobs; obs++) {
      for (int cell = 0; cell < ncell; cell++) {
        int idx = cell + obs * ncell;
        double p = p_step[idx];
        double p_bg = p + bg_rate - p * bg_rate;
        p_step[idx] = p_bg;          // overwrite in place — original value no longer needed
        col_totals[obs] += p_bg;
      }
    }

    // Second pass: normalize and write into the next time-slice.
    for (int obs = 0; obs < nobs; obs++) {
      double total = col_totals[obs];
      if (total <= 0.0) total = 1.0;  // guard against division by zero for degenerate cases
      for (int cell = 0; cell < ncell; cell++) {
        int idx = cell + obs * ncell;
        current[offset_out + idx] = p_step[idx] / total;
      }
    }
  }

  // Attach dimensions so Rcpp returns as 3D array, matching c(ncell, nobs, nsteps)
  current.attr("dim") = IntegerVector::create(ncell, nobs, nsteps);
  return current;
}

// Compute full path propagation negative log-likelihood in one shot.
// Replaces old R chain env_function > apply_kernel > propagate_cpp > log_likelihood
// Parameters:
// par[1:n_env] - environmental attraction parameters (betas)
// par[n_env+1] - sigmoid intercept for environmental attraction
// par[npar-1]  - log(k_exp), exponential movement kernel parameter
// par[npar]    - logit(bg_rate), background movement rate

// [[Rcpp::export]]
double path_propagation_ll_cpp(
  NumericVector par,             // length of npar
  NumericMatrix env_i,
  IntegerMatrix nbhd_i,
  IntegerVector to_dest_vec,
  IntegerVector obs,
  IntegerVector outliers,
  NumericVector inner_dists,
  int ncell_local,
  int n_obs,
  int n_steps,
  int npar,
  int n_env) {
    
    int n_total = ncell_local * n_obs;
    int ncol_inner = nbhd_i.ncol();
    int center = ncell_local / 2;

    double k_exp = std::exp(par[npar - 2]);
    double bg_rate = 1.0 / (1.0 + std::exp(-par[npar - 1]));

    // 1. Environmental attraction (replaces env_function)
    std::vector<double> attract_raw(n_total);
    for (int k = 0; k < n_total; k++) {
      double linear = par[n_env];
      for (int j = 0; j < n_env; j++) {
        linear += par[j] * env_i(k, j);
      }
      attract_raw[k] = 1.0 / (1.0 + std::exp(linear));
    }

    // 2. Dispersal kernel (replaces calculate_dispersal_kernel)
    std::vector<double> kernel(ncol_inner);
    for (int j = 0; j < ncol_inner; j++) {
      kernel[j] = k_exp * std::exp(-k_exp * inner_dists[j]);
    }

    // 3. Transition probabilities (replaces env_function + apply_kernel)
    std::vector<double> attract(n_total * ncol_inner, 0.0);
    for (int row = 0; row < n_total; row++) {
      double row_sum = 0.0;
      for (int col = 0; col < ncol_inner; col++) {
        int idx = nbhd_i(row, col); // 1-indexed from R
        if (IntegerVector::is_na(idx)) continue; // skip NA neighbors
        double val = attract_raw[idx - 1] * kernel[col]; // convert to 0-index
        attract[row + col * n_total] = val; 
        row_sum += val;
      }
      if (row_sum > 0.0) {
        double inv = 1.0 / row_sum;
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] *= inv; // normalize
        }
      }
    }

    // 4. Propagation and likelihood (replaces propagate_cpp + log_likelihood)
    std::vector<double> current(n_total, 0.0);
    std::vector<double> next_buf(n_total, 0.0);
    for (int i = 0; i < n_obs; i++) {
      current[center + i * ncell_local] = 1.0; // initialize step 0
    }
    std::vector<bool> is_outlier(n_obs, false);
    for (int k = 0; k < outliers.size(); k++) {
      int idx = outliers[k]; 
      if (idx >= 1 && idx <= n_obs) {
        is_outlier[idx - 1] = true; // mark outliers (convert to 0-index)
      }
    }
    double floor_val = 2.220446e-16;  // .Machine$double.eps in R, to avoid log(0)
    std::vector<double> ll_per_step(n_steps, 0.0);

    for (int i = 0; i < n_obs - 1; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue; // skip outliers and missing obs
      int cell = obs[i] - 1; // 0-index
      double p = current[cell + i * ncell_local];
      ll_per_step[0] += std::log(std::max(p, floor_val));
    }

    for (int step = 0; step < n_steps - 1; step++) {
      std::fill(next_buf.begin(), next_buf.end(), 0.0);
      for (int k = 0; k < n_total; k++) {
        double incoming = 0.0;
        for (int j = 0; j < ncol_inner; j++) {
          int v = to_dest_vec[k + j * n_total]; // 1-index
          if (IntegerVector::is_na(v)) continue; 
          int src_flat = v - 1;  // 0-index
          int src_row = src_flat % n_total; 
          int src_col = src_flat / n_total;
          incoming += current[src_row] * attract[src_row + src_col * n_total];
        }
        next_buf[k] = incoming + bg_rate - incoming * bg_rate; // apply background rate
      }

      for (int i = 0; i < n_obs; i++) {
        double col_sum = 0.0;
        int base = i * ncell_local;
        for (int cell = 0; cell < ncell_local; cell++) {
          col_sum += next_buf[base + cell];
        }
        if (col_sum > 0.0) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) {
            next_buf[base + cell] *= inv; // normalize
          }
        }
      }

      std::swap(current, next_buf); // move to next step

      for (int i = 0; i < n_obs - 1; i++) {
        if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue; // skip outliers and missing obs
        int cell = obs[i] - 1; // 0-index
        double p = current[cell + i * ncell_local];
        ll_per_step[step + 1] += std::log(std::max(p, floor_val));
      }
    }

    double best_ll = *std::max_element(ll_per_step.begin(), ll_per_step.end());
    return -best_ll; // return negative log-likelihood for minimization
  }
 