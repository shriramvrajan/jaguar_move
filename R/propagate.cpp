#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Compute full path propagation negative log-likelihood in one shot.
// Replaces R chain env_function > apply_kernel > propagate_cpp > log_likelihood
// Parameters:
// par[1:n_env] - environmental attraction parameters (betas)
// par[n_env+1] - sigmoid intercept for environmental attraction
// par[npar-1]  - log(k_exp), exponential movement kernel parameter
// par[npar]    - logit(bg_rate), background movement rate

// Rcpp:export before each function makes it available to R
// [[Rcpp::export]]
double path_propagation_ll_cpp(
  NumericVector attract_raw,
  double k_exp,
  double bg_rate,
  IntegerMatrix nbhd_i,
  IntegerVector to_dest_vec,
  IntegerVector obs,
  IntegerVector outliers,
  NumericVector inner_dists,
  int ncell_local,
  int n_obs,
  int n_steps) {
    
    int n_total = ncell_local * n_obs;
    int ncol_inner = nbhd_i.ncol();
    int center = ncell_local / 2;

    // Dimension checks 
    if (attract_raw.size() != n_total)
      Rcpp::stop("attract_raw has %d elements, expected %d", 
        attract_raw.size(), n_total);
    if (nbhd_i.nrow() != n_total)
      Rcpp::stop("nbhd_i has %d rows, expected %d", nbhd_i.nrow(), n_total);
    if (to_dest_vec.size() != n_total * ncol_inner)
      Rcpp::stop("to_dest_vec has %d elements, expected %d",
                to_dest_vec.size(), n_total * ncol_inner);
    if (obs.size() != n_obs - 1)
      Rcpp::stop("obs has %d elements, expected %d", obs.size(), n_obs - 1);

    // 1. Dispersal kernel (replaces calculate_dispersal_kernel)
    std::vector<double> kernel(ncol_inner);
    for (int j = 0; j < ncol_inner; j++) {
      kernel[j] = k_exp * std::exp(-k_exp * inner_dists[j]);
    }

    // 2. Transition probabilities (replaces env_function + apply_kernel)
    std::vector<double> attract(n_total * ncol_inner, 0.0);
    for (int row = 0; row < n_total; row++) {
      double row_sum = 0.0;
      for (int col = 0; col < ncol_inner; col++) {
        int idx = nbhd_i(row, col); // 1-indexed from R
        if (IntegerVector::is_na(idx)) continue; // skip NA neighbors
        if (idx < 1 || idx > n_total) continue;  // skip out of bounds indices
        double val = attract_raw[idx - 1] * kernel[col]; // convert to 0-index
        if (!std::isfinite(val)) continue;      // non-finite cell -> 0
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
          if (v < 1 || v > n_total * ncol_inner) continue;  // skip OOB obs
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
 
// Function to run diagnostics for an individual, given model objects. 

// [[Rcpp::export]]
List path_propagation_diagnose_cpp(
  NumericVector attract_raw,
  double k_exp,
  double bg_rate,
  IntegerMatrix nbhd_i,
  IntegerVector to_dest_vec,
  IntegerVector obs,
  IntegerVector outliers,
  NumericVector inner_dists,
  int ncell_local,
  int n_obs,
  int n_steps) {

    int n_total = ncell_local * n_obs;
    int ncol_inner = nbhd_i.ncol();
    int center = ncell_local / 2;

    // Dimension checks 
    if (attract_raw.size() != n_total)
      Rcpp::stop("attract_raw has %d elements, expected %d", 
        attract_raw.size(), n_total);
    if (nbhd_i.nrow() != n_total)
      Rcpp::stop("nbhd_i has %d rows, expected %d", nbhd_i.nrow(), n_total);
    if (to_dest_vec.size() != n_total * ncol_inner)
      Rcpp::stop("to_dest_vec has %d elements, expected %d",
                to_dest_vec.size(), n_total * ncol_inner);
    if (obs.size() != n_obs - 1)
      Rcpp::stop("obs has %d elements, expected %d", obs.size(), n_obs - 1);
      
    // 1. Dispersal kernel
    std::vector<double> kernel(ncol_inner);
    for (int j = 0; j < ncol_inner; j++) {
      kernel[j] = k_exp * std::exp(-k_exp * inner_dists[j]);
    }

    // 2. Transition probabilities
    std::vector<double> attract(n_total * ncol_inner, 0.0);
    for (int row = 0; row < n_total; row++) {
      double row_sum = 0.0;
      for (int col = 0; col < ncol_inner; col++) {
        int idx = nbhd_i(row, col);
        if (IntegerVector::is_na(idx)) continue;
        if (idx < 1 || idx > n_total) continue;
        double val = attract_raw[idx - 1] * kernel[col];
        if (!std::isfinite(val)) continue;  // non-finite cells turn to zeroes
        attract[row + col * n_total] = val;
        row_sum += val;
      }
      if (row_sum > 0.0) {
        double inv = 1.0 / row_sum;
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] *= inv;
        }
      }
    }

    // 3. Propagation; same as ll version but keep the surface at each step
    std::vector<double> current(n_total, 0.0);
    std::vector<double> next_buf(n_total, 0.0);
    for (int i = 0; i < n_obs; i++) {
      current[center + i * ncell_local] = 1.0;
    }

    std::vector<bool> is_outlier(n_obs, false);
    for (int k = 0; k < outliers.size(); k++) {
      int idx = outliers[k];
      if (idx >= 1 && idx <= n_obs) is_outlier[idx - 1] = true;
    }

    double floor_val = 2.220446e-16;
    NumericVector ll_by_step(n_steps);
    // Track per-obs ll at each step: n_steps * (n_obs-1), col-major
    int n_transitions = n_obs - 1;
    std::vector<double> ll_obs_all(n_steps * n_transitions, 0.0);

    // Step 0
    for (int i = 0; i < n_transitions; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
      double p = current[obs[i] - 1 + i * ncell_local];
      double ll_i = std::log(std::max(p, floor_val));
      ll_by_step[0] += ll_i;
      ll_obs_all[0 + i * n_steps] = ll_i;
    }

    // Steps {1..n_steps-1}
    for (int step = 0; step < n_steps - 1; step++) {
      std::fill(next_buf.begin(), next_buf.end(), 0.0);
      for (int k = 0; k < n_total; k++) {
        double incoming = 0.0;
        for (int j = 0; j < ncol_inner; j++) {
          int v = to_dest_vec[k + j * n_total];
          if (IntegerVector::is_na(v)) continue;
          if (v < 1 || v > n_total * ncol_inner) continue;
          int src_flat = v - 1;
          int src_row = src_flat % n_total;
          int src_col = src_flat / n_total;
          incoming += current[src_row] * attract[src_row + src_col * n_total];
        }
        next_buf[k] = incoming + bg_rate - incoming * bg_rate;
      }

      for (int i = 0; i < n_obs; i++) {
        double col_sum = 0.0;
        int base = i * ncell_local;
        for (int cell = 0; cell < ncell_local; cell++) col_sum += next_buf[base + cell];
        if (col_sum > 0.0) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) next_buf[base + cell] *= inv;
        }
      }

      std::swap(current, next_buf);

      for (int i = 0; i < n_transitions; i++) {
        if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
        double p = current[obs[i] - 1 + i * ncell_local];
        double ll_i = std::log(std::max(p, floor_val));
        ll_by_step[step + 1] += ll_i;
        ll_obs_all[(step + 1) + i * n_steps] = ll_i;
      }
    }
    
    int best_step = 0;
    double best_ll = ll_by_step[0];
    for (int s = 1; s < n_steps; s++) {
      if (ll_by_step[s] > best_ll) { best_ll = ll_by_step[s]; best_step = s; }
    }

    // Extract per-obs quantities at best step
    std::vector<double> ll_vals;
    std::vector<int> keep;
    for (int i = 0; i < n_transitions; i++) {
      if (is_outlier[i]) continue;
      keep.push_back(i + 1);
      ll_vals.push_back(IntegerVector::is_na(obs[i]) ?
                          NA_REAL : ll_obs_all[best_step + i * n_steps]);
    }
    int n_keep = (int) keep.size();
    NumericVector ll_obs(ll_vals.begin(), ll_vals.end());
    ll_obs.attr("names") = wrap(keep);

    // Snapshot the full surface at best step — but we need to re-run to that step.
    // Re-propagate to best_step:
    std::fill(current.begin(), current.end(), 0.0);
    for (int i = 0; i < n_obs; i++) {
      current[center + i * ncell_local] = 1.0;
    }
    for (int step = 0; step < best_step; step++) {
      std::fill(next_buf.begin(), next_buf.end(), 0.0);
      for (int k = 0; k < n_total; k++) {
        double incoming = 0.0;
        for (int j = 0; j < ncol_inner; j++) {
          int v = to_dest_vec[k + j * n_total];
          if (IntegerVector::is_na(v)) continue;
          if (v < 1 || v > n_total * ncol_inner) continue;
          int src_flat = v - 1;
          int src_row = src_flat % n_total;
          int src_col = src_flat / n_total;
          incoming += current[src_row] * attract[src_row + src_col * n_total];
        }
        next_buf[k] = incoming + bg_rate - incoming * bg_rate;
      }
      for (int i = 0; i < n_obs; i++) {
        double col_sum = 0.0;
        int base = i * ncell_local;
        for (int cell = 0; cell < ncell_local; cell++) col_sum += next_buf[base + cell];
        if (col_sum > 0.0) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) next_buf[base + cell] *= inv;
        }
      }
      std::swap(current, next_buf);
    }

    // Now "current" holds the surface at best_step; subset to kept columns
    NumericMatrix p_surface(ncell_local, n_keep);
    for (int c = 0; c < n_keep; c++) {
      int i = keep[c] - 1;
      for (int cell = 0; cell < ncell_local; cell++)
        p_surface(cell, c) = current[cell + i * ncell_local];
    }
    p_surface.attr("dimnames") = List::create(R_NilValue, wrap(keep));

    return List::create(
      Named("ll_total")   = -best_ll,
      Named("ll_by_step") = ll_by_step,
      Named("best_step")  = best_step + 1,  // back to R 1-indexing
      Named("ll_obs")     = ll_obs,
      Named("p_surface")  = p_surface
    );
}