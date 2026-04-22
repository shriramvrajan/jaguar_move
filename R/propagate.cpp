#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

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

    // Dimension guards
    if (env_i.nrow() != n_total)
      Rcpp::stop("env_i has %d rows, expected %d", env_i.nrow(), n_total);
    if (env_i.ncol() < n_env)
      Rcpp::stop("env_i has %d cols, expected >= %d", env_i.ncol(), n_env);
    if (nbhd_i.nrow() != n_total)
      Rcpp::stop("nbhd_i has %d rows, expected %d", nbhd_i.nrow(), n_total);
    if (to_dest_vec.size() != n_total * ncol_inner)
      Rcpp::stop("to_dest_vec has %d elements, expected %d",
                to_dest_vec.size(), n_total * ncol_inner);
    if (obs.size() != n_obs - 1)
      Rcpp::stop("obs has %d elements, expected %d", obs.size(), n_obs - 1);

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
 