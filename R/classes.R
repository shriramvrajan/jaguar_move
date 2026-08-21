# Attempt to refactor into object-oriented idiom
library(R6)

## Models ======================================================================

# Distance decay kernel model - no habitat selection, no path dependence
distance_model <- R6Class("distance_model",
  public = list(
    prepare_objects = function(trajectory, max_dist, rdf) {
          nbhd <- make_nbhd(rdf = rdf, i = trajectory, sz = max_dist)
          obs <- sapply(2:nrow(nbhd), function(r) {
            nextcell <- which(nbhd[r - 1, ] == trajectory[r])
            if (length(nextcell) == 0) NA else nextcell
          })
          list(nbhd = nbhd, obs = obs, max_dist = max_dist, outliers = integer(0))
        },

    log_likelihood = function(log_k, objects) {
      kernel <- calculate_dispersal_kernel(
        max_dispersal_dist = objects$max_dist,
        kfun = function(x) dexp(x, exp(log_k)))
      nbhd <- objects$nbhd
      kmat <- matrix(rep(kernel, each = nrow(nbhd)), nrow = nrow(nbhd))
      kmat[is.na(nbhd)] <- 0
      p <- kmat / rowSums(kmat)
      valid <- setdiff(seq_along(objects$obs), objects$outliers)
      like <- vapply(valid, function(i) p[i, objects$obs[i]], numeric(1))
      -sum(log(pmax(like, .Machine$double.eps)), na.rm = TRUE)
    },

    fit = function(trajectory, max_dist, rdf, outliers = integer(0)) {
      objects <- self$prepare_objects(trajectory, max_dist, rdf)
      objects$outliers <- outliers
      opt <- optim(log(1), self$log_likelihood, objects = objects,
                  method = "Brent", lower = -5, upper = 5)
      list(log_k = opt$par, k = exp(opt$par), ll = opt$value,
              convergence = opt$convergence, objects = objects)
        }
    )
  )

# Step selection model class 
# Distance kernel model with habitat selection added, but no path dependence.
#  prepare_objects() - prepare data structures for fitting
#  build_kernel()    - build dispersal kernel and calculate attraction
#  log_likelihood()  - calculate log-likelihood of observed steps
#  fit()             - fit the model to a trajectory
#  diagnose()        - more informative output for debugging
step_selection_model <- R6Class("step_selection_model",
  public = list(  
    ## prepare_objects ---------------------------------------------------------
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE) {
      message("Preparing step selection objects...")
      # Neighborhoods for step selection
      nbhd_ss <- make_nbhd(rdf = rdf, i = trajectory, sz = max_dist)

      # Build observations for step selection
      obs <- sapply(2:nrow(nbhd_ss), function(r) {
        local <- nbhd_ss[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })

      if (sim) {
        unique_env <- na.exclude(unique(as.vector(nbhd_ss)))
        env_table <- scale(rdf[unique_env, 1, drop = FALSE])
        env_ss <- as.vector(env_table)
        names(env_ss) <- row.names(env_table)
      } else {
        unique_env <- na.exclude(unique(as.vector(nbhd_ss)))
        env_ss <- scale(rdf[unique_env, ])
      }
      if (any(is.na(env_ss))) env_ss[which(is.na(env_ss))] <- 0

      nbhd_ss <- matrix(as.character(nbhd_ss), 
          nrow = nrow(nbhd_ss), ncol = ncol(nbhd_ss)) # still necessary?

      message("Step selection objects prepared.")
      list(
        env_ss = env_ss,
        nbhd_ss = nbhd_ss, 
        obs = obs,
        max_dist = max_dist,
        outliers = integer(0)  # set externally if needed
      )
    },

    ## build_kernel ------------------------------------------------------------
    build_kernel = function(par, objects, sim, env_type) {
      # Extract kernel parameters
      k_exp   <- ifelse(sim, exp(par[length(par)]), exp(par[length(par) - 1]))
      bg_rate <- ifelse(sim, 0, plogis(par[length(par)]))
      kernel <- if (sim) {
        calculate_dispersal_kernel(max_dispersal_dist = objects$max_dist,
                                   kfun = function(x) dexp(x, k_exp))
      } else {
        stay_kernel(objects$max_dist, k_exp, plogis(par[length(par) - 2]))
      }
      # Calculate environmental attraction
      attract0 <- env_function(objects$env_ss, par, objects$nbhd_ss, sim = sim,
                               type = env_type)
      attract <- apply_kernel(attract0, kernel, bg_rate = bg_rate)
      return(attract)      
    },

    ## log_likelihood ----------------------------------------------------------
    log_likelihood = function(par, objects, sim, env_type = "1o", debug = TRUE) {
      obs      <- objects$obs
      max_dist <- objects$max_dist
      outliers <- objects$outliers
      
      attract <- self$build_kernel(par, objects, sim, env_type = env_type)
      
      indices <- if (length(outliers) > 0) {
        setdiff(seq_along(obs), outliers)
      } else {
        seq_along(obs)
      }
      like <- sapply(indices, function(i) {
        return(attract[i, obs[i]])
      })
      out <- -sum(log(pmax(like, .Machine$double.eps)), na.rm = TRUE)
      # if (is.infinite(out) || is.na(out)) out <- 0
      return(out)
    },

    ## fit ---------------------------------------------------------------------
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, 
                  outliers = integer(0), env_type = "1o", gao = FALSE,
                  gao_control = list()) {
      objects <- self$prepare_objects(trajectory, max_dist, rdf, sim)
      objects$outliers <- outliers
      # Fit model
      # Starting parameters & bounds defined separately for both models, fix when possible
      lbound <- if (sim) {
        c(rep(-Inf, length(par_start)))
      } else {
        c(rep(-Inf, length(par_start) - 3), -10, -5, -30)
      }

      ubound <- if (sim) {
        c(rep(Inf, length(par_start)))
      } else {
        c(rep(Inf, length(par_start) - 3), 10, 5, -2)
      }

      tryCatch({
        if (gao) {
          par_out <- fit_with_gao(le_func = self$log_likelihood, 
            par_start = par_start, lbound = lbound, ubound = ubound, 
            gao_control = gao_control, objects = objects, sim = sim, 
            env_type = env_type)
        } else {
          par_out <- optim(par_start, self$log_likelihood, objects = objects, 
                          sim = sim,  method = "L-BFGS-B", 
                          env_type = env_type,
                          lower = lbound,
                          upper = ubound,
                        control = list(maxit = 400, factr = 1e9))
                        # Max iterations and tolerance factor
        }
        
        ll <- self$log_likelihood(par_out$par, objects, sim, env_type = env_type)
        return(list(
          par = par_out$par,
          ll = ll,
          convergence = par_out$convergence
          # objects = objects
        ))
      }, error = function(e) {
        message(paste("Step selection fitting error:", e$message))
        return(NA)
      })
    },

    ## diagnose ----------------------------------------------------------------
    diagnose = function(par, objects, env_type, sim = FALSE) {
      attract  <- self$build_kernel(par, objects, sim = sim, env_type = env_type)
      obs      <- objects$obs
      outliers <- objects$outliers
      valid    <- setdiff(seq_along(obs), outliers)

      p_obs <- sapply(valid, function(i) attract[i, obs[i]])
      ll_obs <- log(pmax(p_obs, .Machine$double.eps))

      return(list(
        ll_total = -sum(ll_obs, na.rm = TRUE),
        ll_obs = setNames(ll_obs, valid),
        p_obs  = setNames(p_obs, valid),
        p_surface = attract[valid, , drop = FALSE] # rows = obs
      ))
    }
  )
)

# Path propagation model class
#   prepare_objects() - prepare data structures for fitting
#   build_kernel()    - build the dispersal kernel 
#   log_likelihood()  - calculate log-likelihood of observed path
#   fit()             - fit the model to a trajectory
#   diagnose()
#   dispersal_from()  - simulate dispersal from a starting point
path_propagation_model <- R6Class("path_propagation_model",
  public = list(
    propagation_steps = NULL,  # Should be set externally

    ## prepare_objects ---------------------------------------------------------
    prepare_objects = function(trajectory, max_dist, step_size, rdf, sim = FALSE) {
      message("Preparing path propagation objects...")
      # Extended neighborhoods
      nbhd_0 <- make_nbhd(rdf = rdf, i = trajectory, sz = max_dist)

      # Get observations as cell numbers within full neighborhood
      obs <- sapply(2:nrow(nbhd_0), function(r) {
        local <- nbhd_0[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })
      
      # Build inner neighborhood structure. Each entry in the list is the 
      # immediate neighborhood of each cell in the extended neighborhood,
      # as represented by a cell number of raster env_ras
      all_cells <- unique(as.vector(nbhd_0)) %>% na.exclude %>% as.numeric
      nbhd_local <- make_nbhd(rdf = rdf, i = all_cells, sz = step_size)
      cell_to_local <- setNames(seq_along(all_cells), as.character(all_cells))
      nbhd_list <- lapply(seq_len(nrow(nbhd_0)), function(i) {                
        row_inds <- seq_len(ncol(nbhd_0)) + (i - 1) * ncol(nbhd_0)
        names(row_inds) <- nbhd_0[i, ]
        cells_i <- nbhd_0[i, ]
        local_rows <- cell_to_local[as.character(cells_i)]
        out <- matrix(row_inds[as.character(nbhd_local[local_rows, ])], 
                      nrow = length(row_inds), ncol = ncol(nbhd_local))
        return(out)
      })
      nbhd_i <- do.call(rbind, nbhd_list)
      # Reindexing to link row numbers from nbhd_i to cell numbers in env_ras
      nbhd_0_mat <- nbhd_0
      nbhd_0 <- as.vector(t(nbhd_0))
      message("  Inner neighborhood structure built.")
      
      nbhd_flat <- as.integer(nbhd_i)
      valid <- which(!is.na(nbhd_flat) & 
        nbhd_flat >= 1L & nbhd_flat <= nrow(nbhd_i)) # valid entries in nbhd_i
      dest_ids <- nbhd_flat[valid]

      # Sort positions by destination — this groups them
      ord <- order(dest_ids)
      dest_sorted <- dest_ids[ord]
      pos_sorted  <- valid[ord]

      # Run-length encode to get group boundaries
      rle_d <- rle(dest_sorted)
      ncol_td <- ncol(nbhd_i)

      # Cap each group to ncol_td entries, build column indices
      col_idx <- unlist(lapply(pmin(rle_d$lengths, ncol_td), seq_len))

      # Corresponding row indices (repeat each dest id by its capped length)
      capped_lengths <- pmin(rle_d$lengths, ncol_td)
      row_idx <- rep(rle_d$values, capped_lengths)

      # Trim pos_sorted to match: drop excess entries per group
      keep <- unlist(lapply(rle_d$lengths, function(n) {
        c(rep(TRUE, min(n, ncol_td)), rep(FALSE, max(0, n - ncol_td)))
      }))

      to_dest <- matrix(NA_integer_, nrow = nrow(nbhd_i), ncol = ncol_td)
      to_dest[cbind(row_idx, col_idx)] <- pos_sorted[keep]
      dest <- matrix(0, nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
      message("  to_dest and dest matrices built.")

      ## Scale environmental variables 
      unique_cells <- unique(nbhd_0) %>% na.exclude %>% as.numeric
      if (sim) {
        env_table <- scale(rdf$sim1[unique_cells])
        env_s <- as.vector(env_table)
        names(env_s) <- as.character(unique_cells)
        env_i <- env_s[as.character(nbhd_0)]
      } else {
        env_s <- scale(rdf[unique_cells, ])
        env_i <- env_s[match(nbhd_0, unique_cells), ]
      }
      if (any(is.na(env_i))) env_i[which(is.na(env_i))] <- 0

      # Precompute Euclidean distance for each cell in inner neighborhood
      offsets_row <- rep(-step_size:step_size, each = 2 * step_size + 1)
      offsets_col <- rep(step_size:-step_size, times = 2 * step_size + 1)
      inner_dists <- sqrt(offsets_row^2 + offsets_col^2)

      # Return path propagation specific objects
      result <- list(
        env_i = env_i,
        nbhd_i = nbhd_i,
        to_dest = to_dest,
        to_dest_vec = as.vector(to_dest), 
        dest = dest,
        obs = obs,
        max_dist = max_dist,
        step_size = step_size,
        outliers = integer(0),  # set externally if needed
        inner_dists = inner_dists,
        nbhd_0 = nbhd_0_mat
        # mu_env = attributes(env_i)[[2]],
        # sd_env = attributes(env_i)[[3]]
      )
      message("Path propagation objects prepared.")   
      return(result)
    },

    ## log_likelihood ----------------------------------------------------------
    log_likelihood = function(par, objects, sim, debug = FALSE, env_type = "1o") {
      n_obs       <- length(objects$obs) + 1
      ncell_local <- (2 * objects$max_dist + 1)^2
      p_stay  <- if (sim) NA else plogis(par[length(par) - 2])
      k_exp   <- if (sim) exp(par[length(par)]) else exp(par[length(par) - 1])
      bg_rate <- if (sim) 0 else plogis(par[length(par)])
      attract_raw <- as.numeric(env_function(objects$env_i, par, nbhd = NULL,
                                              sim = sim, type = env_type))
      path_propagation_ll_cpp(
        p_stay = p_stay, k_exp = k_exp, bg_rate = bg_rate, 
        attract_raw = attract_raw, nbhd_i = objects$nbhd_i, 
        to_dest_vec = as.integer(objects$to_dest_vec),
        obs = as.integer(objects$obs), outliers = as.integer(objects$outliers),
        multipliers = as.integer(objects$multipliers),
        inner_dists = objects$inner_dists, ncell_local = as.integer(ncell_local),
        n_obs = as.integer(n_obs), n_steps = as.integer(self$propagation_steps))
    },

    ## fit ---------------------------------------------------------------------
    fit = function(trajectory, max_dist, step_size, rdf, par_start, sim = FALSE, 
                  outliers = integer(0), env_type = "1o", multipliers = NULL,
                  gao = FALSE, gao_control = list()) {
      if (is.null(multipliers)) multipliers <- rep(1, length(trajectory) - 1)
      multipliers <- as.integer(multipliers)
      valid_m <- setdiff(seq_along(multipliers), outliers)
      max_m   <- if (length(valid_m) > 0) max(multipliers[valid_m]) else 1

      if (step_size * max_m > max_dist) {
        message(paste0("step_size ", step_size, " x max multiplier ", max_m,
                       " exceeds max_dist ", max_dist, "; skipping."))
        return(NULL)
      }
      depth <- max(1, ceiling(max_dist / step_size))
      cap <- if (sim) depth + 1 else 9
      self$propagation_steps <- min(depth + 1, cap)

      objects <- self$prepare_objects(trajectory, max_dist, step_size, rdf, sim)
      objects$outliers    <- outliers
      objects$multipliers <- multipliers

      trace_log <- list()
      trace_ll <- function(par, objects, sim, lbound, ubound, env_type) {
        ll <- self$log_likelihood(par, objects, sim, env_type = env_type)
        trace_log[[length(trace_log) + 1]] <<- c(par, ll = ll)
        return(ll)
      }

      # Bounds for parameters
      lbound <- if (sim) {
        c(rep(-Inf, length(par_start)))
      } else {
        c(rep(-Inf, length(par_start) - 3), -10, -5, -30)
      }
      ubound <- if (sim) {
        c(rep(Inf, length(par_start)))
      } else {
        c(rep(Inf, length(par_start) - 3), 10, 5, -2)
      }
      # Fit model
      tryCatch({
        par_out <- if (gao) {
          fit_with_gao(self$log_likelihood, par_start, lbound, ubound,
                       gao_control, objects = objects, sim = sim,
                       env_type = env_type)
        } else {
          optim(par_start, trace_ll, objects = objects,                         
                sim = sim,  env_type = env_type, method = "L-BFGS-B", 
                lower = lbound,
                upper = ubound,
                control = list(maxit = 400, factr = 1e9))
        }
                        
        ll <- self$log_likelihood(par_out$par, objects, sim, env_type = env_type)
        traced_ll <- do.call(rbind, trace_log) %>% as.data.frame

        return(list(
          par         = par_out$par,
          ll          = ll,
          convergence = par_out$convergence,
          traced_ll   = traced_ll
          # objects = objects_
        ))
        
      }, error = function(e) {
        message(paste("Path propagation fitting error:", e$message))
        return(NA)
      })
    },
    
    ## diagnose ----------------------------------------------------------------
    diagnose = function(par, objects, sim = FALSE, env_type) {
      ncell_local <- (2 * objects$max_dist + 1)^2
      n_obs <- length(objects$obs) + 1
      p_stay  <- if (sim) NA_real_ else plogis(par[length(par) - 2])
      k_exp <- exp(par[length(par) - 1])
      bg_rate <- plogis(par[length(par)])
      attract_raw <- as.numeric(env_function(objects$env_i, par, nbhd = NULL,
                                sim = sim, type = env_type))

      # debug check
      if (all(is.finite(attract_raw)) & is.finite(k_exp) & is.finite(bg_rate)) {
        print("debug check ✅")
      } else {
        stop()
      }
      path_propagation_diagnose_cpp(
        p_stay = p_stay, k_exp = k_exp, bg_rate = bg_rate,
        attract_raw = attract_raw, nbhd_i = objects$nbhd_i,
        to_dest_vec = as.integer(objects$to_dest_vec),
        obs = as.integer(objects$obs), outliers = as.integer(objects$outliers),
        multipliers = as.integer(objects$multipliers),
        inner_dists = objects$inner_dists, 
        ncell_local = as.integer(ncell_local), n_obs = as.integer(n_obs),
        n_steps = as.integer(self$propagation_steps)
      )
    },

    ## dispersal_from ----------------------------------------------------------
    dispersal_from = function(init_point, par, step_size, n_steps, rdf = brdf,
                              threshold = 1e-6, env_type) {
      
      # FIX BG RATE ISSUE
      k_exp <- exp(par[length(par)]) %>% as.numeric()
      max_displacement <- step_size * n_steps

      # Make it sparse?
      p_current <- setNames(1.0, as.character(init_point))

      for (i in 1:n_steps) {
        message(paste("Step:", i, "of", n_steps))
        p_next <- numeric(0)
        
        for (cell in names(p_current)) {
          p <- p_current[cell]
          if (p < threshold) next

          cell <- as.numeric(cell)
          nbhd <- make_nbhd(rdf = rdf, i = cell, sz = step_size)[1, ]
          nbhd <- nbhd[!is.na(nbhd)]
          
          dists <- sqrt((rdf$row[nbhd] - rdf$row[cell])^2 + 
                        (rdf$col[nbhd] - rdf$col[cell])^2)
          dist_kernel <- dexp(dists, k_exp) 
          env_kernel  <- env_function(rdf[nbhd, 1:6], par, sim = FALSE, type = env_type)
          weights     <- dist_kernel * env_kernel /
                           sum(dist_kernel * env_kernel, na.rm = TRUE)
          
          for (i in seq_along(nbhd)) {
            name         <- as.character(nbhd[i])
            if (!(name %in% names(p_next))) p_next[name] <- 0
            p_next[name] <- p_next[name] + p * weights[i]
          }
        }

        p_current <- p_next
        p_current <- p_current[p_current >= threshold]
      }

      nbhd_full <- make_nbhd(rdf = rdf, i = init_point, sz = max_displacement)
      p_vector <- numeric(length(nbhd_full))
      for (name in names(p_current)) {
        cell <- as.numeric(name)
        idx  <- which(nbhd_full == cell)
        if (length(idx) > 0) p_vector[idx] <- p_current[name]
      }
      p_vector <- p_vector / sum(p_vector, na.rm = TRUE)

      return(list(probs = p_vector,
                  nbhd = nbhd_full,
                  par = par,
                  n_steps = n_steps,
                  step_size = step_size
      ))
    },

    ## stationary_surface ------------------------------------------------------
    stationary_surface = function(par, region_cells, rdf = brdf, step_size = 1,
                                  scale_from = NULL, env_type) {
      k_exp <- exp(par[length(par) - 1]) %>% as.numeric()
      bg_rate <- plogis(par[length(par)]) %>% as.numeric()

      # Scale environment exactly like in prepare_objects
      ref <- if (is.null(scale_from)) region_cells else scale_from
      env_scaled <- scale(rdf[region_cells, 1:6],
                      center = colMeans(rdf[ref, 1:6]),
                      scale  = apply(rdf[ref, 1:6], 2, sd))
      env_scaled[is.na(env_scaled)] <- 0

      # Precompute phi for region cells
      phi_all <- env_function(env_scaled, par, nbhd = NULL, sim = FALSE, 
                              type = env_type)
      names(phi_all) <- as.character(region_cells)

      cell_to_index <- setNames(seq_along(region_cells), as.character(region_cells))
      region_set <- as.character(region_cells)

      ti <- c()  # row indices (from)
      tj <- c()  # column indices (to)
      tx <- c()  # transition probs
      for (k in seq_along(region_cells)) {
        message(paste("Calculating transitions for cell", k, "of", length(region_cells)))
        cell <- region_cells[k]
        nbhd <- make_nbhd(rdf = rdf, i = cell, sz = step_size)[1, ]
        nbhd <- nbhd[!is.na(nbhd)]

        in_region <- as.character(nbhd) %in% region_set
        nbhd <- nbhd[in_region]
        if (length(nbhd) == 0) next

        # same weights as dispersal_from
        dists <- sqrt((rdf$row[nbhd] - rdf$row[cell])^2 + 
                      (rdf$col[nbhd] - rdf$col[cell])^2)
        dist_kernel <- dexp(dists, k_exp)
        env_kernel  <- phi_all[as.character(nbhd)]
        weights     <- dist_kernel * env_kernel /
                       sum(dist_kernel * env_kernel, na.rm = TRUE)
        if (any(is.na(weights))) weights[is.na(weights)] <- 0
        sum_w <- sum(weights)
        if (sum_w == 0) next
        weights <- weights / sum_w
        weights <- weights + bg_rate - weights * bg_rate
        weights <- weights / sum(weights, na.rm = TRUE)

        ti <- c(ti, rep(k, length(nbhd)))
        tj <- c(tj, cell_to_index[as.character(nbhd)])
        tx <- c(tx, weights)
      }

      N <- length(region_cells)
      T_mat <- Matrix::sparseMatrix(i = ti, j = tj, x = tx, dims = c(N, N))
      eig <- tryCatch({
        RSpectra::eigs(t(T_mat), k = 1, which = "LM")
      }, error = function(e) {
        return(NULL) # Return NULL if no convergence
      })
      if (is.null(eig) || length(eig$values) < 1) {
        message("Skipping individual: Eigenvalues failed to converge.")
        return(NA)
      }
      pi_vec <- abs(Re(eig$vectors[, 1]))
      pi_vec <- pi_vec / sum(pi_vec, na.rm = TRUE)
      setNames(pi_vec, as.character(region_cells))
      return(pi_vec)
    }
  )
)

## Theoretical simulations =====================================================

# Configuration parameters for theoretical simulations
#   $max_dist() - maximum dispersal distance for neighborhood
#   $propagation_steps() - number of steps to propagate for path propagation
#   $save() - save configuration parameters to file
simulation_config <- R6Class("simulation_config",
    public <- list(
        # Landscape parameters
        env_size = NULL,
        autocorr_strength = NULL,
        autocorr_range = NULL,
        b_density = NULL,
        b_width = NULL,
        b_length = NULL,
        b_value = NULL,

        # Model parameters
        env_response = NULL,
        step_size    = NULL,
        gen_step     = NULL,
        obs_interval = NULL,

        # Simulation parameters
        n_steps = NULL,
        n_individuals = NULL,
        n_cores = NULL,

        # Output parameters
        name = NULL, 

        initialize = function(name = "default", autocorr_range = 50, 
                        n_individuals = 10, env_size = 200, 
                        autocorr_strength = 5, n_cores = 15,
                        env_response = c(-1.5, 1.5, -0.2, exp(1)),
                        b_density = 0, b_width = 2, b_length = 20, b_value = 99,
                        step_size = 1, gen_step = NULL, obs_interval = 1, 
                        n_steps = 1000) {
            self$name <- name
            self$env_size <- env_size
            self$autocorr_range <- autocorr_range
            self$autocorr_strength <- autocorr_strength
            self$env_response <- env_response
            self$step_size <- step_size
            self$gen_step <- gen_step %||% step_size
            self$obs_interval <- obs_interval
            self$n_individuals <- n_individuals
            self$n_steps <- n_steps
            self$n_cores <- n_cores
            self$b_density <- b_density
            self$b_length <- b_length
            self$b_width <- b_width
            self$b_value <- b_value
        },

        # Derived properties
        max_dist = function() self$gen_step * (self$obs_interval + 1),
        propagation_steps = function() {
          ceiling(self$max_dist() / self$step_size) + 1
        },

        save = function(filepath) {
            params <- list(
                name              = self$name,
                env_size          = self$env_size,
                autocorr_range    = self$autocorr_range,
                autocorr_strength = self$autocorr_strength,
                n_individuals     = self$n_individuals,
                n_steps           = self$n_steps,
                env_response      = self$env_response,
                step_size         = self$step_size,
                gen_step          = self$gen_step,
                obs_interval      = self$obs_interval
            )
            saveRDS(params, filepath)            
        }
    ))

# Movement simulator for path generation and model fitting
#   $simulate_paths() - simulate movement paths based on config
#   $fit_models()    - fit both models to the simulated paths 
movement_simulator <- R6Class("movement_simulator",
    public = list(
        config = NULL,

        initialize = function(config) {
            self$config <- config
        },

        ## simulate_paths ------------------------------------------------------
        simulate_paths = function() {  # terra not working well with parallel
            message(paste("Simulating", self$config$n_individuals, 
                            "paths for", self$config$name))
            
            # Extract config values
            env_size <- self$config$env_size
            autocorr_strength <- self$config$autocorr_strength
            autocorr_range <- self$config$autocorr_range
            env_response <- self$config$env_response
            n_individuals <- self$config$n_individuals
            n_steps <- self$config$n_steps
            step_size <- self$config$step_size
        
            results <- list()
            for (i in 1:self$config$n_individuals) {
              config_seed <- sum(utf8ToInt(self$config$name))
              set.seed(i + 4001 + config_seed)
              land_i <- gen_landscape(
                  size = self$config$env_size,
                  s = self$config$autocorr_strength,
                  r = self$config$autocorr_range,
                  b_density = self$config$b_density,
                  b_length = self$config$b_length,
                  b_width = self$config$b_width,
                  b_value = self$config$b_value
              )
              rdf <- raster_to_df(land_i)
              message(paste0("Path #: ", i, " / ", self$config$n_individuals))
              x0 <- ceiling(self$config$env_size / 2)
              y0 <- ceiling(self$config$env_size / 2)
              
              path_i <- jag_path(x0, y0, self$config$n_steps, 
                                par = self$config$env_response, 
                                neighb = self$config$gen_step,
                                rdf = rdf)

              # Proportion of barriers encountered
              nbhd <- make_nbhd(rdf = rdf, i = path_i$cell, sz = self$config$gen_step)
              barrier_cells <- which(rdf$sim1 == self$config$b_value)
              prop_barrier <- sum(apply(nbhd, 1, function(x) sum(x %in% barrier_cells))) / length(nbhd)

              results[[i]] <- list(path = path_i, rdf = rdf, 
                                   connect = calc_connectivity(land_i,
                                 barrier_threshold = 0.9 * self$config$b_value),
                                   prop_barrier = prop_barrier)
            }            
            return(results)
    },
    
    ## fit_models --------------------------------------------------------------
    fit_models = function(paths, parallel = TRUE, n_cores = NULL) {
        message(paste("Fitting models for", self$config$name))
        
        # Extract config information
        n_cores <- n_cores %||% self$config$n_cores
        max_dist <- self$config$max_dist()
        propagation_steps <- self$config$propagation_steps()
        par_start <- c(-1.5, 1.5, -0.2, log(1.0))

        obs_interval <- self$config$obs_interval
        step_size <- self$config$step_size
        n_individuals <- self$config$n_individuals
        env_size <- self$config$env_size
        autocorr_strength <- self$config$autocorr_strength
        autocorr_range <- self$config$autocorr_range
        b_density <- self$config$b_density
        b_length <- self$config$b_length
        b_width <- self$config$b_width
        b_value <- self$config$b_value

        ss_model <- step_selection_model$new()
        pp_model <- path_propagation_model$new()
        
        if (parallel) {
          registerDoParallel(n_cores) 
          # Get ALL functions from global environment
          all_functions <- ls(globalenv())[sapply(ls(globalenv()), 
                                              function(x) is.function(get(x)))]
          
          results <- foreach(i = seq_len(n_individuals),
            .packages = c("gstat"),
            .export = c(  # Model classes
                          "step_selection_model", "path_propagation_model",
                          "ss_model", "pp_model",
                          # All functions
                          all_functions,
                          # Data
                          "paths",
                          # Parameters
                          "max_dist", "propagation_steps", "par_start", 
                          "step_size", "env_size", "autocorr_strength", 
                          "autocorr_range", "obs_interval", "b_density", 
                          "b_length", "b_width", "b_value")) %dopar% {            

            message(paste0("Fitting individual #: ", i))

            path_i <- paths[[i]]$path
            rdf    <- paths[[i]]$rdf
            trajectory_df <- cbind(path_i$x, path_i$y)
            ind <- seq(1, nrow(trajectory_df), obs_interval + 1)
            trajectory_df <- trajectory_df[ind, ]
            nrow_ras <- max(rdf$row)
            ncol_ras <- max(rdf$col)
            trajectory <- (trajectory_df[, 1] - 1) * ncol_ras + trajectory_df[, 2]
            # trajectory <- (nrow_ras - trajectory_df[, 2]) * ncol_ras + trajectory_df[, 1]

            # Build full neighborhood
            nbhd_full <<- make_nbhd(rdf = rdf, i = seq_len(nrow(rdf)), 
                              sz = step_size)
            
            # Get environmental mean and sd for local neighborhood
            visited_cells <- unique(as.vector(make_nbhd(rdf = rdf, 
                                              i = trajectory, sz = max_dist)))
            visited_cells <- visited_cells[!is.na(visited_cells)]
            env_mu <- mean(rdf[visited_cells, 1])
            env_sd <- sd(rdf[visited_cells, 1])

            # Fit both models
            ss_result <- ss_model$fit(trajectory, max_dist, rdf = rdf, 
                        par_start, sim = TRUE)
            
            pp_model$propagation_steps <- propagation_steps
            pp_result <- pp_model$fit(trajectory, max_dist, rdf = rdf,
              par_start, sim = TRUE, step_size = step_size)
            
            return(list(
              step_selection = ss_result,
              path_propagation = pp_result,
              path = path_i,
              rdf = rdf,
              env_mu = env_mu,
              env_sd = env_sd
            ))
          }
        } else { 
          results <- list()
          for (i in 1:self$config$n_individuals) {
              message(paste0("Fitting individual #: ", i))
              
              path_i <- paths[[i]]$path
              rdf    <- paths[[i]]$rdf
              # Prepare trajectory
              trajectory_df <- cbind(path_i$x, path_i$y)
              ind <- seq(1, nrow(trajectory_df), self$config$obs_interval + 1)
              trajectory_df <- trajectory_df[ind, ]
              nrow_ras <- max(rdf$row)
              ncol_ras <- max(rdf$col)
              trajectory <- (trajectory_df[, 1] - 1) * ncol_ras + trajectory_df[, 2]
              # trajectory <- (nrow_ras - trajectory_df[, 2]) * ncol_ras + trajectory_df[, 1]

              # Build individual neighborhood
              nbhd_full <<- make_nbhd(rdf = rdf, i = seq_len(nrow(rdf)), 
                                    sz = self$config$step_size)
              
              # Get environmental mean and sd for local neighborhood
              visited_cells <- unique(as.vector(make_nbhd(rdf = rdf, 
                                                i = trajectory, sz = max_dist)))
              visited_cells <- visited_cells[!is.na(visited_cells)]
              env_mu <- mean(rdf[visited_cells, 1])
              env_sd <- sd(rdf[visited_cells, 1])

              # Fit both models
              ss_result <- ss_model$fit(trajectory, max_dist, rdf = rdf, 
                         par_start = par_start, sim = TRUE)
              pp_result <- pp_model$fit(trajectory, max_dist, rdf = rdf, 
                         par_start = par_start, sim = TRUE, step_size = step_size)
              results[[i]] <- list(
                step_selection = ss_result,
                path_propagation = pp_result,
                path = path_i,
                rdf = rdf,
                env_mu = env_mu,
                env_sd = env_sd
              )
          }
        }
      
      return(results)
    },

    # plot_fits ----------------------------------------------------------------
    plot_fits = function(results, save = FALSE) {
      for (i in seq_along(results)) {
        print(i)
        result <- results[[i]]
        # print(result)
        if (!is.list(result) || identical(result$step_selection, NA) ||
              identical(result$path_propagation, NA)) next

        if (save) cairo_pdf(filename = paste0("figs/sims/currentpaths/fit_", 
          self$config$name, "_ind_", i, ".pdf"), width = 10, height = 5)
        par(mfrow = c(1, 2))
        
        # Panel 1: path on landscape
        path_obs <- result$path[seq(1, nrow(result$path), 
                                self$config$obs_interval + 1), ]
        land_mat <- matrix(result$rdf[[1]], 
                   nrow = self$config$env_size, 
                   ncol = self$config$env_size, byrow = TRUE)
        land_ras <- rast(land_mat)
        ext(land_ras) <- c(0.5, self$config$env_size + 0.5,
                          0.5, self$config$env_size + 0.5)
        buffer <- 10
        xlim <- c(max(1, range(result$path$x)[1] - buffer), 
                  min(self$config$env_size, range(result$path$x)[2] + buffer))
        ylim <- c(max(1, range(result$path$y)[1] - buffer), 
                  min(self$config$env_size, range(result$path$y)[2] + buffer))
        image(1:self$config$env_size, 1:self$config$env_size, land_mat,
                xlim = xlim, ylim = ylim,
                main = paste0("#", i), col = hcl.colors(50, "viridis"))
        points(path_obs$x, path_obs$y, col = rgb(1, 1, 1, 0.3), 
                pch = 19, cex = 0.5)
        lines(path_obs$x, path_obs$y, col = rgb(1, 1, 1, 0.3), lwd = 1)

        # terra::plot(land_ras, xlim = xlim, ylim = ylim)
        # points(path_obs$x, path_obs$y, col = rgb(1, 1, 1, 0.3), 
        #         pch = 19, cex = 0.5)
        # lines(path_obs$x, path_obs$y, col = rgb(1, 1, 1, 0.3), lwd = 1)
        
        # Panel 2: environmental response curves
        visible <- make_nbhd(rdf = result$rdf, i = unique(result$path$cell),
                            sz = self$config$max_dist()) %>% as.vector()
        visible_env <- result$rdf$sim1[visible]
        if (any(visible_env == self$config$b_value)) {
          visible_env <- visible_env[-which(visible_env == self$config$b_value)]
        }
        env_used <- result$rdf$sim1[result$path$cell]
        
        env_range <- range(visible_env)
        x_vals <- seq(env_range[1], env_range[2], 0.1)
        y_true <- 1 / (1 + exp(self$config$env_response[1] + 
                              self$config$env_response[2] * x_vals + 
                              self$config$env_response[3] * x_vals^2))
        y_true <- y_true / max(y_true)  # Normalize 
        plot(x_vals, y_true, type = "l", col = "black", lwd = 2,
            xlim = env_range,  # Use actual data range
            xlab = "Environmental variable", ylab = "Relative attraction")

        points(visible_env, jitter(rep(0.05, length(visible_env)), amount = 0.02), 
             col = rgb(0.7, 0.7, 0.7, 0.5), pch = 16, cex = 0.5)
        points(env_used, jitter(rep(0.1, length(env_used)), amount = 0.02), 
             col = rgb(0, 0, 0, 0.5), pch = 16, cex = 0.5)

        x_scaled <- (x_vals - result$env_mu) / result$env_sd
        # Fitted curves
        if (!is.na(result$step_selection[1])) {
          y_ss <- 1 / (1 + exp(result$step_selection$par[1] + 
                              result$step_selection$par[2] * x_scaled + 
                              result$step_selection$par[3] * x_scaled^2))
          y_ss <- y_ss / max(y_ss)
          lines(x_vals, y_ss, col = "red", lwd = 2)
        }
        if (!is.na(result$path_propagation[1])) {
          y_pp <- 1 / (1 + exp(result$path_propagation$par[1] + 
                              result$path_propagation$par[2] * x_scaled + 
                              result$path_propagation$par[3] * x_scaled^2))
          y_pp <- y_pp / max(y_pp)
          lines(x_vals, y_pp, col = "blue", lwd = 2)
        }
        abline(v = x_vals[which.max(y_true)], lty = 2, col = "black")
        legend("topright", c("True", "SS", "PP"), col = c("black", "red", "blue"), lwd = 2)

        dev.off()      
      }
    }
))

# Batch runner to coordinate multiple simulations
#   $autocorr_range_study() - create batch for fragmentation study
#   $run_all() - run all simulations in the batch
#   $get_results() - extract performance comparison for plotting
#   $plot_fragmentation_effect() - plot effect on model performance
simulation_batch <- R6Class("simulation_batch",
  public = list(
    configs = list(),
    results = list(),
    done = NA,
    
    # run_all ==================================================================
    run_all = function(parallel = TRUE, n_cores = NULL) {

      self$results <- list()  
      for (config in self$configs) {
        if (config$name %in% self$done) {
          message(paste("Skipping, already run:", config$name))
          self$results[[config$name]] <- readRDS(
            paste0("data/output/sim_", config$name, ".rds")
          )
          next
        }
        message(paste("Running simulation:", config$name))
        # Create simulator
        sim <- movement_simulator$new(config)
        n_cores <- n_cores %||% config$n_cores
        paths <- sim$simulate_paths()                           # Generate paths
        fit_results <- sim$fit_models(paths, parallel, n_cores)     # Fit models
        sim$plot_fits(fit_results)                   # Plot fits for diagnostics
        self$results[[config$name]] <- list(                     # Store results
            config = config,
            paths = paths,
            fits = fit_results
        ) 
        saveRDS(self$results[[config$name]], 
                paste0("data/output/sim_", config$name, ".rds"))
        message(paste("Completed:", config$name))
      }
      return(self$results)
    },
    
    # get_results ==============================================================
    get_results = function() {
      results_list <- list()
      summary_list <- list()
      for (name in names(self$results)) {
        print(paste("Processing results for:", name))
        fits <- self$results[[name]]$fits
        config <- self$results[[name]]$config

        # Extract log-likelihoods for each individual
        n_obs <- floor(config$n_steps / (config$obs_interval + 1))
        ll_step <- sapply(fits, function(x) {
          if (is.na(x[1]) || length(x[1]) == 0) NA else x$step_selection$ll
        })
        ll_pp <- sapply(fits, function(x) {
          if (is.na(x[2]) || length(x[2]) == 0) NA else x$path_propagation$ll
        })
        ll_diff <- sapply(seq_along(ll_step), function(x) {
          if (is.na(ll_step[x]) || is.na(ll_pp[x]) || ll_pp[x] == 0 || ll_step[x] == 0) {
            return(NA) 
          } else {
            return(ll_step[x] - ll_pp[x])
          }
        })
        ll_diff_per_obs <- ll_diff / (n_obs - 1)  # Per observation

        # Connectivity metric
        connectivity <- sapply(seq_along(fits), function(i) {
          paths_i <- self$results[[name]]$paths
          if (!is.null(paths_i[[i]]$connect)) paths_i[[i]]$connect else NA
        })

        # Proportion of barriers encountered
        prop_barrier <- sapply(seq_along(fits), function(i) {
          paths_i <- self$results[[name]]$paths
          if (!is.null(paths_i[[i]]$prop_barrier)) paths_i[[i]]$prop_barrier else NA
        })

        # Calculate AIC difference (positive means path propagation better)
        aic_diff <- 2 * (ll_step - ll_pp) # - 4? (2 parameters, 2k - 2LL)
        
        results_list[[name]] <- data.frame(
          config_name = name,
          individual = seq_along(fits),
          gen_step = config$gen_step,
          step_size = config$step_size,
          ll_step = ll_step,
          ll_pp = ll_pp,
          ll_diff_per_obs = ll_diff_per_obs,
          connect = connectivity,
          prop_barrier = prop_barrier,
          autocorr_range = config$autocorr_range,
          obs_interval = config$obs_interval,
          b_density = config$b_density,
          aic_diff = aic_diff,
          n_successful = sum(!is.na(aic_diff))
        )
      }

      results_df <- do.call(rbind, results_list)
      configs    <- unique(results_df$config_name)
      summary_df <- data.frame(
        config_name = configs,
        autocorr_range = sapply(configs, function(x) {  
            unique(results_df$autocorr_range[results_df$config_name == x])
        }),
        gen_step = sapply(configs, function(x) {
            unique(results_df$gen_step[results_df$config_name == x])
        }),
        obs_interval = sapply(configs, function(x) {
            unique(results_df$obs_interval[results_df$config_name == x])
        }),
        b_density = sapply(configs, function(x) {
            unique(results_df$b_density[results_df$config_name == x])
        }),
        ll_diff_per_obs = sapply(configs, function(x) {
            median(results_df$ll_diff_per_obs[results_df$config_name == x],
                  na.rm = TRUE)
        }))
      return(list(results = results_df, summary = summary_df))
    }
))

## Empirical testing ===========================================================

# Jaguar class to handle data related to a single individual
#   $get_track()     - get processed track data with movement metrics
#   $get_landscape() - get cropped landscape raster and dataframe for the track
jaguar <- R6Class("jaguar",
  public = list(
    id            = NULL,
    obs_interval  = NULL,
    track         = NULL,
    track_cells   = NULL,
    landscape     = NULL,
    results       = NULL,
    multipliers   = NULL,
    outliers      = NULL,
    max_dist      = NULL,
    max_multiple  = NULL,
    scale_time    = NULL,

    initialize = function(id = NULL, results = NULL, max_multiple = 2,
                          obs_interval = 0, scale_time = FALSE) {
      self$id <- as.numeric(id)
      self$max_multiple <- max_multiple
      self$scale_time <- scale_time
      
      self$obs_interval <- obs_interval
      self$track <- self$get_track()
      self$track_cells <- cellFromXY(brazil_ras, 
                          self$track[, c("longitude", "latitude")])
      self$landscape <- self$get_landscape()
      self$results <- results[which(results$ID == self$id), ] %>%
                        as.data.frame

      dt_base   <- get_mode(round(na.exclude(self$track$dt)))
      dt_scaled <- self$track$dt[-1] / dt_base
      dt_discrete <- round(dt_scaled)
      self$outliers <- which(dt_discrete < 1 | dt_discrete > self$max_multiple)
      self$multipliers <- if (self$scale_time) {
                              as.integer(pmax(dt_discrete, 1L)) 
                            } else { 
                              rep(1L, length(dt_discrete))
                            }

      sl_all <- self$track$sl
      sl_emp <- na.exclude(if (length(self$outliers)) 
         sl_all[-(self$outliers + 1)] else sl_all)
      self$max_dist <- ceiling(1.1 * max(sl_emp) / 1000)
    },

    ## get_track ---------------------------------------------------------------
    get_track = function() {
      dat <- jag_move[ID == self$id]
      keep <- seq(1, nrow(dat), by = self$obs_interval + 1)
      dat <- dat[keep, ]
      dat$timestamp <- as.POSIXct(dat$timestamp, 
                              format = "%m/%d/%Y %H:%M")
      dat$year <- as.numeric(format(dat$timestamp, "%Y"))
      dat$year <- ifelse(dat$year > 23, dat$year + 1900, dat$year + 2000)

      dat$mon <- as.numeric(format(dat$timestamp, "%m"))
      dat$day <- as.numeric(format(dat$timestamp, "%d"))
      dat$hr <- format(dat$timestamp, "%H:%M")
      dat$hr <- as.numeric(gsub(":[0-9][0-9]", "", dat$hr))
      
      path <- jag_move[ID == self$id][keep, ]
      path$t <- lubridate::mdy_hm(as.character(path$timestamp))
      path <- vect(path, geom = c("longitude", "latitude"), crs = wgs84)
      path <- project(path, epsg5880)
      path <- track(x = crds(path)[, 1], y = crds(path)[, 2], 
                    t = path$t, id = path$ID, crs = epsg5880)
      st <- steps(path)
      
      dat$sl <- c(NA, st$sl_)             # step lengths in m
      dat$ta <- c(NA, st$ta_)             # turn angles in radians
      dat$dir <- c(NA, st$direction_p)    # bearing in radians
      dat$dt <- c(NA, as.numeric(st$dt_, units = "mins")) # time interval in minutes
      dat$spd <- dat$sl / dat$dt
      return(dat[, c("timestamp", "longitude", "latitude", "ID", "year", "mon", 
                    "day", "hr", "sl", "ta", "dir", "dt", "spd")])
    },

    ## get_landscape -----------------------------------------------------------
    get_landscape = function() {
      # Bounding box of the jaguar's trajectory with 0.1 degree buffer
      track <- self$get_track()
      track_extent <- ext(terra::vect(track, 
                          geom = c("longitude", "latitude"), crs = wgs84))
      track_extent <- ext(
        track_extent[1] - 0.1,
        track_extent[2] + 0.1,
        track_extent[3] - 0.1,
        track_extent[4] + 0.1
      )
      brazil_ras_crop <- terra::crop(brazil_ras, track_extent)
      brazil_ras_crop_df <- raster_to_df(brazil_ras_crop)
      return(list(raster = brazil_ras_crop, dataframe = brazil_ras_crop_df))
    },

    # eval_ss ------------------------------------------------------------------
    eval_ss = function(par, env_type = "1o") {
      ss <- step_selection_model$new()
      obj <- ss$prepare_objects(self$track_cells, self$max_dist, rdf = brdf)
      obj$outliers <- self$outliers
      ss$diagnose(par, obj, env_type = env_type)
    },

    # eval_pp ------------------------------------------------------------------
    # Use case 2: evaluate LL at given parameters
    eval_pp = function(par, n_jump, env_type = "1o") {
      if (n_jump == "max") n_jump = self$max_dist - 1
      setup <- self$pp_objects(n_jump, env_type)
      if (is.null(setup)) stop("n_jump too coarse to resolve the multipliers.")
      ll  <- setup$model$log_likelihood(par, setup$objects, sim = FALSE, env_type = env_type)
      return(setup$model$diagnose(par, setup$objects, env_type = env_type))
    },

    # fit_ss --------------------------------------------------------------
    fit_ss = function(par_start, env_type = "1o", clamp = Inf, gao = FALSE,
                      gao_control = list()) {
      env_idx <- seq_len(length(par_start) - 3)
      par_start[env_idx] <- pmin(pmax(par_start[env_idx], -clamp), clamp)

      ss  <- step_selection_model$new()
      res <- ss$fit(self$track_cells, self$max_dist, rdf = brdf,
                    par_start = par_start, outliers = self$outliers,
                    env_type = env_type, gao = gao, gao_control = gao_control)
      if (length(res) == 1 && is.na(res[1])) return(NULL)
      res
    },
    
    # fit_pp -------------------------------------------------------------------
    # Use case 1: fit PP with n_jump
    fit_pp = function(par_start, env_type = "1o", max_jump = NULL, 
               n_jump = NULL, clamp = 2, gao = FALSE, gao_control = list()) {
      max_dist <- self$max_dist
      if (is.null(max_jump)) max_jump <- min(max_dist - 1, 8)
      print(max_jump)

      env_idx <- seq_len(length(par_start) - 3)
      par_start[env_idx] <- pmin(pmax(par_start[env_idx], -clamp), clamp)

      best <- list(result = NULL, ll = Inf, n_jump = 0)
      worse_count <- 0

      for (n_j in 0:max_jump) {
        if (worse_count > 2) break
        if (!is.null(n_jump) && n_jump != n_j) next
        message(paste0("n_jump = ", n_j))
        pp <- path_propagation_model$new()
        res <- pp$fit(self$track_cells, max_dist, rdf = brdf,
          par_start = par_start, outliers = self$outliers,
          step_size = n_j + 1, env_type = env_type,
          multipliers = self$multipliers, gao = FALSE, gao_control = gao_control)
        if (is.null(res) || (length(res) == 1 && is.na(res[1]))) next
        par_start <- res$par                             # warm-start next n_j
        if (res$ll < best$ll) {
          best <- list(result = res, ll = res$ll, n_jump = n_j)
          worse_count <- 0
        } else {
          worse_count <- worse_count + 1
        }
      }

      if (!is.null(best$result) && gao) {
        nj <- best$n_jump
        message(paste0("Refining n_jump = ", nj, " with GAO"))
        pp <- path_propagation_model$new()
        pp$propagation_steps <- ceiling(max_dist / (nj + 1))
        ref <- pp$fit(self$track_cells, max_dist, rdf = brdf,
          par_start = best$result$par, outliers = self$outliers,
          step_size = nj + 1, env_type = env_type,
          multipliers = self$multipliers,
          gao = TRUE, gao_control = gao_control)
        if (!is.null(ref) && !is.na(ref[1]) && ref$ll < best$ll) {
          best <- list(result = ref, ll = ref$ll, n_jump = nj)
        }
      }
      best$result$n_jump <- best$n_jump
      return(best$result)
    },

    # pp_objects ---------------------------------------------------------------
    pp_objects = function(n_jump, env_type = "1o", holdout_frac = NULL, training = TRUE) {
      tc   <- self$track_cells
      out  <- self$outliers
      mult <- self$multipliers
      max_dist <- self$max_dist

      if (!is.null(holdout_frac)) {
        hold <- seq_len(ceiling(length(tc) * holdout_frac))
        if (training) {
          tc   <- tc[hold]
          mult <- mult[seq_len(length(hold) - 1)]
        } else {
          tc   <- tc[-hold]
          out  <- out[out > max(hold)] - length(hold)
          mult <- mult[(length(hold) + 1):length(mult)]
        }
      }

      step_size <- n_jump + 1
      pp <- path_propagation_model$new()
      valid <- which(!(seq_along(mult) %in% out))
      max_m <- if (length(valid) > 0) max(mult[valid]) else 1L
      if (step_size * max_m > max_dist) return(NULL)
      depth <- max(1, ceiling(max_dist / step_size))
      pp$propagation_steps <- min(depth + 1, 9)
      obj <- pp$prepare_objects(tc, max_dist, step_size, rdf = brdf)
      obj$outliers    <- as.integer(out)
      obj$multipliers <- as.integer(mult)
      list(model = pp, objects = obj, max_dist = max_dist)
    },

    # holdout_pp ---------------------------------------------------------------
    # Use case 3: holdout analysis
    holdout_pp = function(par_start, env_type = "1o", frac = 0.7, max_jump = NULL,
                          clamp = 2) {
      max_dist <- self$max_dist
      if (is.null(max_jump)) max_jump <- min(max_dist - 1, 8)

      env_idx <- seq_len(length(par_start) - 3)
      par_start[env_idx] <- pmin(pmax(par_start[env_idx], -clamp), clamp)

      # Fit on training set
      tc_train   <- self$track_cells[seq_len(ceiling(length(self$track_cells) * frac))]
      best <- list(result = NULL, ll = Inf, n_jump = 0)
      worse_count <- 0

      for (n_j in 0:max_jump) {
        if (worse_count > 2) break
        setup <- self$pp_objects(n_j, env_type, holdout_frac = frac, training = TRUE)
        res <- setup$model$fit(tc_train, max_dist, rdf = brdf,
          par_start = par_start, outliers = setup$objects$outliers,
          step_size = n_j + 1, env_type = env_type,
          multipliers = as.integer(setup$objects$multipliers))
        if (!is.na(res[1])) par_start <- res$par
        if (!is.na(res[1]) && res$ll < best$ll) {
          best <- list(result = res, ll = res$ll, n_jump = n_j)
          worse_count <- 0
        } else {
          worse_count <- worse_count + 1
        }
      }

      # Evaluate on holdout set
      hold_setup <- self$pp_objects(best$n_jump, env_type, holdout_frac = frac,
                                    training = FALSE)
      ll_hold <- hold_setup$model$log_likelihood(best$result$par, hold_setup$objects,
        sim = FALSE, env_type = env_type)

      return(list(train = best$result, n_jump = best$n_jump,
          holdout_nll = ll_hold, holdout_n = length(hold_setup$objects$obs)))
    },

    # calculate_ll ------------------------------------------------------------
    calculate_ll = function(env_type, par_ss = NULL, par_pp = NULL, n_jump = NULL) {
      if (is.null(self$results)) stop("No fitted results.")
      npar <- length(grep("ss_par", names(self$results)))
      ss_ind <- grep("ss_par", names(self$results))
      pp_ind <- grep("pp_par", names(self$results))
      if (is.null(par_ss)) par_ss <- as.numeric(self$results[, ss_ind])
      if (is.null(par_pp)) par_pp <- as.numeric(self$results[, pp_ind])
      if (is.null(n_jump))  n_jump <- ifelse(is.na(self$results$pp_njump), 0, 
                                              self$results$pp_njump)

      return(list(ss = self$eval_ss(par_ss, env_type = env_type),
          pp = self$eval_pp(par_pp, n_jump = n_jump, env_type = env_type)))
    },

    # calculate_null_ll -------------------------------------------------------
    calculate_null_ll = function() {
      sl     <- self$track$sl
      ocheck <- length(self$outliers) >= 1
      max_dist <- self$max_dist

      fit <- distance_model$new()$fit(
        self$track_cells, max_dist, rdf = brdf, outliers = self$outliers)

      obj   <- fit$objects
      valid <- setdiff(seq_along(obj$obs), obj$outliers)
      valid <- valid[!is.na(obj$obs[valid])]
      m_i   <- rowSums(!is.na(obj$nbhd))[valid]

      return(list(ll_unif = sum(log(m_i)),
                  ll_dist = fit$ll,
                  k       = fit$k,
                  n       = length(valid)))
    },

    # calculate_limit_case_ll -------------------------------------------------
    calculate_limit_case_ll = function() {
      ocheck <- length(self$outliers) >= 1
      max_dist <- self$max_dist
      ## ... incomplete
    },

    # benchmark ---------------------------------------------------------------
    benchmark = function(n_jump = 0, env_type = "1o") {
      track <- self$track
      track_cells <- self$track_cells
      outliers    <- self$outliers
      multipliers <- self$multipliers
      max_dist    <- self$max_dist
      
      step_size <- n_jump + 1
      pp <- path_propagation_model$new()
      pp$propagation_steps <- min(max(1, ceiling(max_dist / step_size)) + 1, 9)
      objects <- pp$prepare_objects(track_cells, max_dist, step_size, rdf = brdf)
      objects$outliers <- outliers
      objects$multipliers <- as.integer(multipliers)
      
      gc()
      mem_mb <- gc()[2, 2]
      cat("Memory after prepare:", mem_mb, "MB\n")
      cat("osize:", jag_meta$osize[jag_meta$ID == self$id], "MB\n")
      cat("Multiplier:", mem_mb / (jag_meta$osize[jag_meta$ID == self$id]), "\n")
      cat("ncell_local:", (2 * max_dist + 1)^2, " n_obs:", length(track_cells), "\n")
      
      return(list(pp = pp, objects = objects, max_dist = max_dist, outliers = outliers,
          multipliers = multipliers))
    }
  ))

# Batch runner for empirical data fitting
# Needs a bit of refactoring for cyclomatic complexity
#   $run_all() - run fitting for all individuals in config
#   $process_individual() - process a single individual 
empirical_batch <- R6Class("empirical_batch",
  public = list(
    config = NULL,
    k_fit = NULL,
    ss_warm_par = NULL,
    out_dir = NULL,
    results = list(),

    initialize = function(config, k_fit, ss_warm_par) {
      self$config <- config
      self$k_fit  <- k_fit
      self$ss_warm_par <- ss_warm_par
      self$out_dir <- file.path("data/output", sprintf("%s_%s_m%s_o%s%s",
        switch(config$model_type, "ss", "pp"), config$env_type,
        config$max_multiple %||% 1, config$obs_interval %||% 0,
        if (config$holdout_set) "_holdout" else ""))
      dir.create(self$out_dir, recursive = TRUE, showWarnings = FALSE)
    },

    ## run_all -----------------------------------------------------------------
    run_all = function() {
      # Determine which individuals to process
      if (is.null(self$config$individuals)) {
        individuals_to_process <- jag_id$jag_id
      } else {
        individuals_to_process <- self$config$individuals
      }
      
      # Check for already completed results, then order remainder in reverse
      # order of total object size (proxy for time it takes to run) so that 
      # we don't waste core time
      i_todo <- setdiff(individuals_to_process, private$get_completed_ids())
      i_todo <- i_todo[order(jag_meta$osize[match(i_todo, jag_meta$ID)])]
      # smallest first
      message(paste0("Processing ", length(i_todo), " individuals"))
      
      core_cap <- self$config$n_cores

      waves <- if (self$config$model_type == 1) {
        # SS: negligible per-worker memory, no batching
        list(i_todo)
      } else {
        # PP wave optimization
        # cost proxy per individual: total object size
        mem_budget <- self$config$mem_budget
        per_core_extra <- 1.5e9  # ~1.5gb overhead per core?
        inflation_risk <- 2      # how much does osize inflate, conservative estimate
        effective_budget <- mem_budget - core_cap * per_core_extra
        cost <- pmin(2e9, jag_meta$osize[match(i_todo, jag_meta$ID)] * inflation_risk)

        out <- list()
        w <- integer(0)
        wcost <- 0

        for (j in seq_along(i_todo)) {
          if (cost[j] > effective_budget) {
            message(paste0("Individual ", i_todo[j], " exceeds budget alone (",
                           signif(cost[j] / 1e9, 3), " GB); will run solo."))
          } else if (length(w) > 0 && (wcost + cost[j] > effective_budget ||
                                       length(w) >= core_cap)) {
            out[[length(out) + 1]] <- w
            w <- integer(0)
            wcost <- 0
          }
          w <- c(w, i_todo[j])
          wcost <- wcost + cost[j]
        }
        if (length(w) > 0) out[[length(out) + 1]] <- w
        out
      }

      message(paste0("Scheduled ", length(i_todo), " individuals in ", 
                      length(waves), " waves"))

      if (self$config$parallel) {
        batch_config <- self$config
        k_fit <- self$k_fit
        ss_warm <- self$ss_warm_par

        cl <- makeCluster(core_cap, outfile = "")
        on.exit(try(stopCluster(cl), silent = TRUE), add = TRUE)
        registerDoParallel(cl)
        parallel::clusterEvalQ(cl, {
          source("R/functions.R")
          source("R/classes.R")
        })

        for (w_i in seq_along(waves)) {
          wave <- waves[[w_i]]
          message(paste0("Wave ", w_i, "/", length(waves), ": ",
                    length(wave), " individuals on", core_cap, " cores"))

          foreach(i = wave, .export = c("batch_config", "k_fit", "ss_warm"),
                  .errorhandling = "pass") %dopar% {
                    ss_model <- step_selection_model$new()
                    pp_model <- path_propagation_model$new()
                    batch <- empirical_batch$new(batch_config, 
                                               k_fit = k_fit, ss_warm = ss_warm)
                    batch$process_individual(i, ss_model, pp_model)
                  } 
          message("Parallel processing complete")    
        }
        stopCluster(cl)
      } else {
          results <- list()
          for (i in i_todo) {                    
            ss_model <- step_selection_model$new()
            pp_model <- path_propagation_model$new()
            results[[length(results) + 1]] <- self$process_individual(i, ss_model, pp_model)
          }
      }

      # If there are already completed results, load and combine
      results <- lapply(seq_len(82), function(i) {
        f <- file.path(self$out_dir, paste0("out_", jag_id$jag_id[i], ".rds"))
        if (file.exists(f)) return(readRDS(f)) else return(NA)
      })
      self$results <- results
      return(results)
    },
    
    ## process_individual ------------------------------------------------------
    process_individual = function(i, ss_model, pp_model, diagnose_par = NULL) {
      message(paste0("Processing jaguar #", i))
      # maybe just skip jaguar 114 or something idk 

      # Create jaguar instance and get data
      mm <- self$config$max_multiple %||% 1
      jag <- jaguar$new(i, max_multiple = mm, 
                        scale_time   = self$config$scale_time %||% FALSE,
                        obs_interval = self$config$obs_interval %||% 0)
      track_cells <- jag$track_cells
      outliers    <- jag$outliers
      multipliers <- jag$multipliers
      max_dist    <- jag$max_dist
            
      # Holdout set mode
      if (self$config$holdout_set && nrow(jag$track) > 100) {
        hold <- seq_len(ceiling(nrow(jag$track) * self$config$holdout_frac))
        if (self$config$fit_model) {
          track_cells <- track_cells[hold]
          multipliers <- multipliers[seq_len(length(hold) - 1)]
        } else {
          track_cells <- track_cells[-hold]
          outliers <- outliers[outliers > max(hold)]
          if (length(outliers) > 0) outliers <- outliers - length(hold)
          multipliers <- multipliers[(length(hold) + 1):length(multipliers)] 
        }
      }

      n_valid <- length(track_cells) - 1 - 
                    sum(outliers %in% seq_len(length(track_cells) - 1))
      if (n_valid < 30) { # minimum sample size, usable transitions
        message(paste0(i, ": sample size ", n_valid, " after thinning; skip"))
        saveRDS(NA, file.path(self$out_dir, paste0("NA_", i, ".rds")))
        return(NA)
      }

      # Diagnosis mode: skip fitting and return diagnostic info
      if (!is.null(diagnose_par)) {
        step_size <- diagnose_par$n_jump + 1
        pp_model$propagation_steps <- min(max(1, ceiling(max_dist / step_size)), 8)
        env_type <- diagnose_par$env_type

        obj_ss <- ss_model$prepare_objects(track_cells, max_dist, rdf = brdf)
        obj_ss$outliers <- outliers
        obj_pp <- pp_model$prepare_objects(track_cells, max_dist, step_size, 
                                          rdf = brdf)
        obj_pp$outliers <- as.integer(outliers)
        obj_pp$multipliers <- as.integer(multipliers)

        print(diagnose_par$pp)

        return(list(
          ss = ss_model$diagnose(diagnose_par$ss, obj_ss, env_type = env_type),
          pp = pp_model$diagnose(diagnose_par$pp, obj_pp, env_type = env_type)
        ))
      }

      # Starting parameters
      ss_i <- if (is.list(self$ss_warm_par)) {
                self$ss_warm_par[[which(jag_id$jag_id == i)]]
              } else NULL
      if (is.list(ss_i) && length(ss_i$par) == self$config$npar) {
        par_start <- ss_i$par
        message("Warm starting from step selection fit")
      } else {
        par_start <- c(rnorm(self$config$npar - 3, sd = 0.1), # environmental vars
                      -0.5,                                      # logit p_stay
                      log(self$k_fit$k[which(self$k_fit$id == i)]),   # k_exp 
                      -20)                                         # bg_rate ~ 0  
      }

      # Fit models based on config
      if (self$config$model_type == 1) {
        # Step selection
        result <- ss_model$fit(track_cells, max_dist, rdf = brdf, 
                               env_type = self$config$env_type,
                               par_start = par_start, outliers = outliers,
                               gao = self$config$gao %||% FALSE,
                              gao_control = self$config$gao_control %||% list())
      } else if (self$config$model_type == 2) {
        check_file <- file.path(self$out_dir, paste0("partial_", i, ".rds"))
        # The check_file loop prevents rerunning n_jump iterations for the case
        # where a core crashes without completing evaluation (happens).

        checkpoint <- if (file.exists(check_file)) {
          message("Resuming from checkpoint.")
          readRDS(check_file)
        } else {
          list(best_result = NULL, best_ll = Inf, best_n_jump = 0, 
               worse_count = 0, done_jumps = integer(0))
        }

        try_jump <- function(n_j) {
          pp_model$propagation_steps <- ceiling(max_dist / (n_j + 1))
          message(paste0("n_jump = ", n_j, ", steps = ", pp_model$propagation_steps))
          res <- pp_model$fit(track_cells, max_dist, rdf = brdf, 
            par_start = par_start, outliers = outliers, step_size = n_j + 1,
            env_type = self$config$env_type, multipliers = multipliers,
            gao = FALSE, # no GAO to find jump, just later for final fit
            gao_control = self$config$gao_control %||% list())
          if (!is.na(res[1])) par_start <<- res$par # warm start next radius
          if (!is.na(res[1]) && res$ll < checkpoint$best_ll) {
            checkpoint$best_ll <<- res$ll
            checkpoint$best_result <<- res
            checkpoint$best_n_jump <<- n_j
            checkpoint$worse_count <<- 0
          } else if (!is.null(checkpoint$best_result)) {
            checkpoint$worse_count <<- checkpoint$worse_count + 1 
          }
          checkpoint$done_jumps <<- c(checkpoint$done_jumps, n_j)
          saveRDS(checkpoint, check_file)
        }

        max_jump <- min(max_dist - 1, 8)
        for (n_jump in 0:max_jump) {
          if (n_jump %in% checkpoint$done_jumps || checkpoint$worse_count > 2) {
            next # If past two have been consecutively worse, skip to end
          } else {
            try_jump(n_jump)
          }
        }
        # if (!((max_dist - 1) %in% checkpoint$done_jumps)) try_jump(max_dist - 1)
        # For 2-consecutive-worse case, also evaluate maximum n_jump. 
        # Necessary for model comparison, since it's the special case 
        # where path propagation is the same as step selection.

        result <- checkpoint$best_result
        if (!is.null(result) && self$config$gao) {
          nj <- checkpoint$best_n_jump
          message(paste0("Running best n_jump = ", nj, " through GAO"))
          pp_model$propagation_steps <- ceiling(max_dist / (nj + 1))
          ref <- pp_model$fit(track_cells, max_dist, rdf = brdf,
            par_start = result$par, outliers = outliers, step_size = nj + 1,
            env_type = self$config$env_type, multipliers = multipliers,
            gao = TRUE, gao_control = self$config$gao_control %||% list())
          if (!is.na(ref[1]) && ref$ll < result$ll) result <- ref
        }
        if (!is.null(result)) result$n_jump <- checkpoint$best_n_jump
        if (!is.null(result)) result$n_jump <- checkpoint$best_n_jump
        if (file.exists(check_file)) file.remove(check_file)
      }
      # save individual result
      if (is.list(result)) result$n_valid <- n_valid # sample size of transitions
      save_ok <- is.list(result) && !is.null(result$par)
      saveRDS(if (save_ok) result else NA,
              file.path(self$out_dir, paste0(if (save_ok) "out_" else "NA_", i, ".rds")))
      
      return(result)
    }
  ),
  
  private = list(
    get_completed_ids = function() {
      outfiles <- list.files(self$out_dir, pattern = "^(out|NA)_.*\\.rds$") 
      ids <- gsub("^(out|NA)_|\\.rds$", "", outfiles) %>% as.numeric
      if (any(is.na(ids))) ids <- ids[!is.na(ids)]
      return(ids)
    }
  ))

results_set <- R6Class("results_set",
  public = list(
    r_ss = NULL,
    r_pp = NULL,
    env_type = NULL,
    res_table = NULL,

    initialize = function(r_ss = NULL, r_pp = NULL, env_type) {
      if (!is.null(r_ss))
        self$r_ss <- if (is.character(r_ss)) readRDS(r_ss) else r_ss
      if (!is.null(r_pp)) {
        self$r_pp <- if (is.character(r_pp)) readRDS(r_pp) else r_pp
        self$r_pp <- lapply(self$r_pp, function(x) {
          if (is.list(x)) { x$traced_ll <- NULL; x } else x
        })
      }
      self$env_type <- env_type
      self$res_table <- self$results_table()
    },

    ## results_table -----------------------------------------------------------
    results_table = function() {
      tabulate_model <- function(results, prefix, has_njump = FALSE) {
        first <- Position(is.list, results)
        if (is.na(first)) return(NULL)
        npar <- length(results[[first]]$par)
        ncols <- npar + 3L + has_njump

        mat <- t(vapply(results, function(r) {
          if (!is.list(r)) return(rep(NA, ncols))
          c(r$par, r$ll, r$convergence, r$n_valid %||% NA,
            if (has_njump) r$n_jump %||% NA)
        }, numeric(ncols)))

        cnames <- c(paste0(prefix, "_par", seq_len(npar)),
                    paste0(prefix, c("_ll", "_conv", "_n")))
        if (has_njump) cnames <- c(cnames, paste0(prefix, "_njump"))

        k <- ifelse(has_njump, npar + 2, npar)
        aic <- ifelse(is.na(mat[, npar + 1]), NA, 2 * mat[, npar + 1] + 2 * k)

        out <- as.data.frame(mat)
        names(out) <- cnames
        out[[paste0(prefix, "_aic")]] <- aic
        # results[[i]] ~ jag_id$jag_id[i], realign to jag_meta rows
        ids <- if (!is.null(names(results))) as.numeric(names(results))
                else as.numeric(jag_id$jag_id)
        out <- out[match(as.numeric(jag_meta$ID), ids), , drop = FALSE]
        rownames(out) <- NULL
        return(out)
      }

      parts <- list(jag_meta[, c("ID", "biome")])
      if (!is.null(self$r_ss))
        parts[[length(parts) + 1]] <- tabulate_model(self$r_ss, "ss")
      if (!is.null(self$r_pp))
        parts[[length(parts) + 1]] <- tabulate_model(self$r_pp, "pp", 
                                                      has_njump = TRUE)
      do.call(cbind, parts)
    }
  ))