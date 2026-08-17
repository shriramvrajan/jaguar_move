#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
#include <limits>
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
  double p_stay,
  IntegerMatrix nbhd_i,
  IntegerVector to_dest_vec,
  IntegerVector obs,
  IntegerVector outliers,
  IntegerVector multipliers,   // time-multiplier per transition (length n_obs-1)
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

    // 1. Dispersal kernel, centre held out at p_stay
    std::vector<double> kernel(ncol_inner);
    int kcenter = 0;
    for (int j = 1; j < ncol_inner; j++)
      if (inner_dists[j] < inner_dists[kcenter]) kcenter = j;

    double off_sum = 0.0;
    for (int j = 0; j < ncol_inner; j++) {
      kernel[j] = (j == kcenter) ? 0.0 : k_exp * std::exp(-k_exp * inner_dists[j]);
      off_sum += kernel[j];
    }
    if (off_sum > 0.0) {
      double scale = (1.0 - p_stay) / off_sum;
      for (int j = 0; j < ncol_inner; j++) kernel[j] *= scale;
    }
    kernel[kcenter] = p_stay;

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
      if (std::isfinite(row_sum) && row_sum >= std::numeric_limits<double>::min()) {
        double inv = 1.0 / row_sum;
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] *= inv;
        }
      } else {
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] = 0.0;
        }
      }
    }

    // 4. Propagation & likelihood (replaces old propagate_cpp + log_likelihood)
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

    // Find max multiplier among valid (non-outlier, non-NA) transitions
    int max_m = 1;
    for (int i = 0; i < n_obs - 1; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
      int m_i = multipliers[i];
      if (m_i > max_m) max_m = m_i;
    }
    int max_r = (n_steps - 1) / max_m;   // max base-rate that stays in range
    std::vector<double> total_ll(max_r + 1, 0.0);  // one accumulator per candidate r

    // Step 0: depth=0 for all transitions, contributes to r=0 (the SS baseline)
    for (int i = 0; i < n_obs - 1; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
      double p = current[obs[i] - 1 + i * ncell_local];
      total_ll[0] += std::log(std::max(p, floor_val));
    }

    // Steps 1..n_steps-1: propagate, then credit each transition to its base-rate
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
        for (int cell = 0; cell < ncell_local; cell++) {
          col_sum += next_buf[base + cell];
        }
        if (std::isfinite(col_sum) && col_sum >= std::numeric_limits<double>::min()) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) next_buf[base + cell] *= inv;
        } else {
          Rcpp::stop("block %d unnormalisable at step %d (col_sum = %g)", i, step, col_sum);
        }
      }

      std::swap(current, next_buf);

      int depth = step + 1;  // current propagation depth after this round
      for (int i = 0; i < n_obs - 1; i++) {
        if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
        int m_i = multipliers[i];
        // This transition contributes to base-rate r iff depth == r * m_i
        if (depth % m_i == 0) {                       // depth is a multiple of m_i
          int r = depth / m_i;                         // the base-rate this depth implies
          if (r >= 1 && r <= max_r) {                  // within searchable range
            double p = current[obs[i] - 1 + i * ncell_local];
            total_ll[r] += std::log(std::max(p, floor_val));
          }
        }
      }
    }

    double best_ll = *std::max_element(total_ll.begin(), total_ll.end());
    return -best_ll;
  }
 
// Function to run diagnostics for an individual, given model objects.
// [[Rcpp::export]]
List path_propagation_diagnose_cpp(
  NumericVector attract_raw,
  double k_exp,
  double bg_rate,
  double p_stay,
  IntegerMatrix nbhd_i,
  IntegerVector to_dest_vec,
  IntegerVector obs,
  IntegerVector outliers,
  IntegerVector multipliers,
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
    int kcenter = 0;
    for (int j = 1; j < ncol_inner; j++)
      if (inner_dists[j] < inner_dists[kcenter]) kcenter = j;   // distance 0
    double off_sum = 0.0;
    for (int j = 0; j < ncol_inner; j++) {
      kernel[j] = (j == kcenter) ? 0.0 : k_exp * std::exp(-k_exp * inner_dists[j]);
      off_sum += kernel[j];
    }
    if (off_sum > 0.0) {
      double scale = (1.0 - p_stay) / off_sum;
      for (int j = 0; j < ncol_inner; j++) kernel[j] *= scale;
    }
    kernel[kcenter] = p_stay;

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
      if (std::isfinite(row_sum) && row_sum >= std::numeric_limits<double>::min()) {
        double inv = 1.0 / row_sum;
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] *= inv;
        }
      } else {
        for (int col = 0; col < ncol_inner; col++) {
          attract[row + col * n_total] = 0.0;
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
        if (std::isfinite(col_sum) && col_sum >= std::numeric_limits<double>::min()) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) next_buf[base + cell] *= inv;
        } else {
          Rcpp::stop("block %d unnormalisable at step %d (col_sum = %g)", i, step, col_sum);
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
    
    // Base-rate profiling (mirrors the fitting function)
    int max_m = 1;
    for (int i = 0; i < n_transitions; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
      if (multipliers[i] > max_m) max_m = multipliers[i];
    }
    int max_r = (n_steps - 1) / max_m;
    NumericVector total_ll_by_r(max_r + 1, 0.0);

    // r=0: all valid transitions at step 0
    for (int i = 0; i < n_transitions; i++) {
      if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
      total_ll_by_r[0] += ll_obs_all[0 + i * n_steps];   // step 0 value
    }
    // r>=1: each transition at step r * m_i
    for (int r = 1; r <= max_r; r++) {
      for (int i = 0; i < n_transitions; i++) {
        if (is_outlier[i] || IntegerVector::is_na(obs[i])) continue;
        int target = r * multipliers[i];
        if (target < n_steps) {
          total_ll_by_r[r] += ll_obs_all[target + i * n_steps];
        }
      }
    }

    int best_r = 0;
    double best_r_ll = total_ll_by_r[0];
    for (int r = 1; r <= max_r; r++) {
      if (total_ll_by_r[r] > best_r_ll) { best_r_ll = total_ll_by_r[r]; best_r = r; }
    }
    
    int best_step = 0;
    double best_ll = ll_by_step[0];
    for (int s = 1; s < n_steps; s++) {
      if (ll_by_step[s] > best_ll) { best_ll = ll_by_step[s]; best_step = s; }
    }

    // Per-transition evaluation depth implied by the chosen base rate
    std::vector<int> eval_depth(n_transitions, 0);
    for (int i = 0; i < n_transitions; i++) eval_depth[i] = best_r * multipliers[i];

    // Per-obs log-likelihoods, read at each transition's own depth
    std::vector<double> ll_vals;
    std::vector<int> keep;
    for (int i = 0; i < n_transitions; i++) {
      if (is_outlier[i]) continue;
      keep.push_back(i + 1);
      int d = eval_depth[i];
      bool bad = IntegerVector::is_na(obs[i]) || d >= n_steps;
      ll_vals.push_back(bad ? NA_REAL : ll_obs_all[d + i * n_steps]);
    }
    int n_keep = (int) keep.size();
    NumericVector ll_obs(ll_vals.begin(), ll_vals.end());
    ll_obs.attr("names") = wrap(keep);

    // Surfaces: re-propagate once, harvesting each column at its own depth
    NumericMatrix p_surface(ncell_local, n_keep);
    int max_depth = 0;
    for (int c = 0; c < n_keep; c++) {
      int d = eval_depth[keep[c] - 1];
      if (d > max_depth) max_depth = d;
    }
    if (max_depth > n_steps - 1) max_depth = n_steps - 1;

    std::fill(current.begin(), current.end(), 0.0);
    for (int i = 0; i < n_obs; i++) current[center + i * ncell_local] = 1.0;

    for (int c = 0; c < n_keep; c++) {          // depth 0 (only if best_r == 0)
      int i = keep[c] - 1;
      if (eval_depth[i] != 0) continue;
      for (int cell = 0; cell < ncell_local; cell++)
        p_surface(cell, c) = current[cell + i * ncell_local];
    }

    for (int step = 0; step < max_depth; step++) {
      std::fill(next_buf.begin(), next_buf.end(), 0.0);
      for (int k = 0; k < n_total; k++) {
        double incoming = 0.0;
        for (int j = 0; j < ncol_inner; j++) {
          int v = to_dest_vec[k + j * n_total];
          if (IntegerVector::is_na(v)) continue;
          if (v < 1 || v > n_total * ncol_inner) continue;
          int src_flat = v - 1;
          int src_row  = src_flat % n_total;
          int src_col  = src_flat / n_total;
          incoming += current[src_row] * attract[src_row + src_col * n_total];
        }
        next_buf[k] = incoming + bg_rate - incoming * bg_rate;
      }
      for (int i = 0; i < n_obs; i++) {
        double col_sum = 0.0;
        int base = i * ncell_local;
        for (int cell = 0; cell < ncell_local; cell++) col_sum += next_buf[base + cell];
        if (std::isfinite(col_sum) && col_sum >= std::numeric_limits<double>::min()) {
          double inv = 1.0 / col_sum;
          for (int cell = 0; cell < ncell_local; cell++) next_buf[base + cell] *= inv;
        } else {
          Rcpp::stop("block %d unnormalisable at step %d (col_sum = %g)", i, step, col_sum);
        }
      }
      std::swap(current, next_buf);

      int depth = step + 1;
      for (int c = 0; c < n_keep; c++) {
        int i = keep[c] - 1;
        if (eval_depth[i] != depth) continue;
        for (int cell = 0; cell < ncell_local; cell++)
          p_surface(cell, c) = current[cell + i * ncell_local];
      }
    }
    p_surface.attr("dimnames") = List::create(R_NilValue, wrap(keep));

    return List::create(
      Named("ll_total")   = -best_r_ll,
      Named("ll_by_step") = ll_by_step,
      Named("ll_by_baserate") = total_ll_by_r, 
      Named("best_step")  = best_step + 1,  // back to R 1-indexing
      Named("best_r")         = best_r, 
      Named("eval_depth") = wrap(eval_depth),
      Named("ll_obs")     = ll_obs,
      Named("p_surface")  = p_surface
    );
}