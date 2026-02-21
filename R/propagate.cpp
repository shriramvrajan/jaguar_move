#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector propagate_cpp(NumericVector attract_vec,  // flattened attract matrix (col-major), length slice_size * ncol_dest
                            IntegerVector to_dest_vec,  // flattened connectivity map (col-major, 1-indexed from R), same length
                            int ncell,             // cells per local neighborhood: (2*max_dist+1)^2
                            int nobs,              // number of observations (trajectory length)
                            int nsteps,            // propagation steps
                            double bg_rate,        // background movement rate (logistic-transformed in R before passing)
                            int ncol_dest) {       // columns in dest/to_dest matrix (= ncol of inner nbhd, typically 9)

  // slice_size is the number of elements in one time-slice of the 3D array.
  // Each slice is a matrix of ncell rows (cells in local neighborhood) by nobs columns (observations).
  // R stores this column-major: element [cell, obs] lives at index cell + obs * ncell.
  int slice_size = ncell * nobs;

  // Full output, flat vector representing 3D array current[cell, obs, step]
  // Stored column-major following R conventions:
  //   current[c, o, s] = current[c + o*ncell + s*slice_size]
  int total_size = slice_size * nsteps;
  NumericVector current(total_size, 0.0);

  // Initialize step 0: probability 1.0 at the center cell of each observation's neighborhood.
  // ncell is always odd ((2k+1)^2), so integer division gives the 0-based center index.
  // This matches R's: center <- ncell_local / 2 + 0.5 (which is 1-indexed).
  int center = ncell / 2;
  for (int obs = 0; obs < nobs; obs++) {
    current[center + obs * ncell] = 1.0;
  }

  // Pre-allocate working memory. Reused every iteration — no allocation inside the loop.
  // p_step[i] will hold the total incoming probability to cell i (across all neighbor directions).
  std::vector<double> p_step(slice_size);

  // One normalization total per observation column.
  std::vector<double> col_totals(nobs);

  // === Main propagation loop: steps 0 through nsteps-2 produce slices 1 through nsteps-1 ===
  for (int j = 0; j < nsteps - 1; j++) {

    // Byte offsets into `current` for the source (j) and destination (j+1) slices.
    int offset_in  = j * slice_size;
    int offset_out = (j + 1) * slice_size;

    // Zero out the accumulator for this iteration.
    std::fill(p_step.begin(), p_step.end(), 0.0);

    // --- Fused gather + rowSum ---
    // In R this was three steps:
    //   step_prob <- as.vector(current[,,j]) * attract[]    (recycle-multiply)
    //   dest[]    <- step_prob[to_dest_vec]                  (gather by index)
    //   p_step    <- rowSums(dest)                           (sum columns)
    //
    // Here we do it in one pass: for each destination cell (row of to_dest),
    // walk across its ncol_dest incoming connections, compute the source
    // probability * transition weight on the fly, and accumulate directly.
    for (int row = 0; row < slice_size; row++) {
      double incoming = 0.0;

      for (int c = 0; c < ncol_dest; c++) {
        // to_dest is stored column-major in R: element [row, c] is at row + c * slice_size.
        int td_val = to_dest_vec[row + c * slice_size];

        // R integer NA comes through as NA_INTEGER (INT_MIN). Skip missing connections.
        if (td_val == NA_INTEGER) continue;

        // Convert from R's 1-indexing to C++ 0-indexing.
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

    // --- Background rate and normalization ---
    // The "noisy OR" combination: p_bg = p + bg - p*bg
    // This ensures every cell retains at least bg_rate probability of being occupied,
    // even if no directed movement reaches it. When bg_rate ≈ 0 (typical: plogis(-15)),
    // this is nearly a no-op but prevents numerical zeros.
    //
    // After mixing, normalize each observation column to sum to 1.

    // First pass: apply bg_rate and compute column totals.
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

  // Attach dimensions so Rcpp returns this as a 3D array, matching R's array(dim = c(ncell, nobs, nsteps)).
  current.attr("dim") = IntegerVector::create(ncell, nobs, nsteps);
  return current;
}