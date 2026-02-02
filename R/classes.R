# Attempt to refactor into object-oriented idiom
library(R6)

## Models ======================================================================

# Step selection model class
#  $prepare_objects() - prepare data structures for fitting
#  $build_kernel()    - build the dispersal kernel and calculate attraction
#  $dispersal_from()  - simulate dispersal from a starting point
#  $log_likelihood()  - calculate log-likelihood of observed steps
#  $fit()             - fit the model to a trajectory
step_selection_model <- R6Class("step_selection_model",
  public = list(  
    ## prepare_objects ---------------------------------------------------------
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE) {
      message("Preparing step selection objects")
      # Neighborhoods for step selection
      nbhd_ss <- make_nbhd(rdf = rdf, i = trajectory, sz = max_dist)

      # Build observations for step selection
      obs <- sapply(2:nrow(nbhd_ss), function(r) {
        local <- nbhd_ss[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })

      if (sim) {
        env_table <- scale(rdf[unique(nbhd_ss), 1, drop = FALSE])
        env_ss <- as.vector(env_table)
        names(env_ss) <- row.names(env_table)
      } else {
        env_ss <- scale(rdf[unique(nbhd_ss), ])
      }
      if (any(is.na(env_ss))) env_ss[which(is.na(env_ss))] <- 0

      nbhd_ss <- matrix(as.character(nbhd_ss), 
          nrow = nrow(nbhd_ss), ncol = ncol(nbhd_ss)) # still necessary?

      list(
        env_ss = env_ss,
        nbhd_ss = nbhd_ss, 
        obs = obs,
        max_dist = max_dist,
        outliers = integer(0)  # set externally if needed
      )
    },

    ## build_kernel ------------------------------------------------------------
    build_kernel = function(par, objects, sim) {
      # Extract kernel parameters
      k_exp   <- ifelse(sim, exp(par[length(par)]), exp(par[length(par) - 1]))
      bg_rate <- ifelse(sim, 0, plogis(par[length(par)]))
      kernel <- calculate_dispersal_kernel(max_dispersal_dist = objects$max_dist, 
                                          kfun = function(x) dexp(x, k_exp))
      # Calculate environmental attraction
      attract0 <- env_function(objects$env_ss, par, objects$nbhd_ss, sim = sim)
      attract <- apply_kernel(attract0, kernel, bg_rate = bg_rate)
      return(attract)      
    },

    ## log_likelihood ----------------------------------------------------------
    log_likelihood = function(par, objects, sim, debug = TRUE) {
      obs      <- objects$obs
      max_dist <- objects$max_dist
      outliers <- objects$outliers
      
      attract <- self$build_kernel(par, objects, sim)
      
      indices <- if (length(outliers) > 0) {
        setdiff(seq_along(obs), outliers)
      } else {
        seq_along(obs)
      }
      like <- sapply(indices, function(i) {
        return(attract[i, obs[i]])
      })
      out <- -sum(log(like), na.rm = TRUE)

      # if (is.infinite(out) || is.na(out)) out <- 0

      print(paste("Current LL:", round(-out, 4)))
      print(paste("SS parameters:", paste(round(par, 4), collapse = ", ")))
      # saveRDS(list(out = out, array = like, par = par, o = objects),
      #         "data/output/current_iter_ss.rds")
      return(out)
    },

    ## fit ---------------------------------------------------------------------
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, 
                  outliers = integer(0)) {
      objects <- self$prepare_objects(trajectory, max_dist, rdf, sim)
      objects$outliers <- outliers
      print(paste("PAR START:", as.character(par_start), ","))
      # Fit model
      # Starting parameters & bounds defined separately for both models, fix when possible
      lbound <- if (sim) {
        c(rep(-Inf, length(par_start)))
      } else {
        c(rep(-Inf, length(par_start) - 2), -5, -30)
      }

      ubound <- if (sim) {
        c(rep(Inf, length(par_start)))
      } else {
        c(rep(Inf, length(par_start) - 2), 5, -2)
      }

      tryCatch({
        par_out <- optim(par_start, self$log_likelihood, objects = objects, 
                        sim = sim,  method = "L-BFGS-B",
                        lower = lbound,
                        upper = ubound,
                        control = list(
                          maxit = 2000,     # More iterations
                          factr = 1e5
                        ))
        
        ll <- self$log_likelihood(par_out$par, objects, sim)
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

    dispersal_from = function(init_point, par, rdf = brdf, raster = brazil_ras,
                              max_dist, sim = FALSE) {
      objects <- self$prepare_objects(init_point, max_dist, rdf, sim, raster)
      attract <- self$build_kernel(par, objects, sim)
      return(list(probs = attract[1, ],
                  nbhd = objects$nbhd_ss,
                  par = par))
    }
  )
)

# Path propagation model class
#   $prepare_objects() - prepare data structures for fitting
#   $build_kernel()    - build the dispersal kernel 
#   $dispersal_from()  - simulate dispersal from a starting point
#   $log_likelihood()  - calculate log-likelihood of observed path
#   $fit()             - fit the model to a trajectory
path_propagation_model <- R6Class("path_propagation_model",
  public = list(
    propagation_steps = 10,  # Default, should be set externally

    ## prepare_objects ---------------------------------------------------------
    prepare_objects = function(trajectory, max_dist, step_size, rdf, sim = FALSE) {
      ## This is where the magic happens. Indexing is everything!
      message("Preparing path propagation objects")   
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
      nbhd_list <- lapply(seq_len(nrow(nbhd_0)), function(i) {                
        row_inds <- seq_len(ncol(nbhd_0)) + (i - 1) * ncol(nbhd_0)
        names(row_inds) <- nbhd_0[i, ]
        out <- matrix(row_inds[as.character(nbhd_full[nbhd_0[i, ], ])], 
                      nrow = length(row_inds), ncol = ncol(nbhd_full))
        return(out)
      })
      nbhd_i <- do.call(rbind, nbhd_list)
      # Reindexing to link row numbers from nbhd_i to cell numbers in env_ras
      nbhd_0 <- as.vector(t(nbhd_0))
      
      # Build connectivity matrix. Rows are each cell of nbhd_i, columns are
      # row numbers from nbhd_i that can be reached in one step. All row lengths
      # standardized to ncol(nbhd_i) with NAs
      to_dest <- tapply(seq_len(length(nbhd_i)), nbhd_i, function(x) {  
        if (length(x) <= ncol(nbhd_i)) {
          out <- c(x, rep(NA, ncol(nbhd_i) - length(x)))
        } else {
          out <- x[seq_len(ncol(nbhd_i))]
        }
        return(out)
      })
      to_dest <- t(matrix(unlist(to_dest), nrow = ncol(nbhd_i), 
                                           ncol = nrow(nbhd_i)))
      dest <- matrix(0, nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
      
      unique_cells <- unique(nbhd_0) %>% na.exclude %>% as.numeric
      # env_s <- if (sim) scale(rdf$sim1[nbhd_0]) else scale(rdf[unique_cells, ])
      # env_i <- env_s[match(nbhd_0, unique_cells), ]
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
        outliers = integer(0)  # set externally if needed
        # mu_env = attributes(env_i)[[2]],
        # sd_env = attributes(env_i)[[3]]
      )
      
      return(result)
    },

    ## build_kernel ------------------------------------------------------------
    build_kernel = function(par, objects, sim) { 
      # Numerical values
      max_dist   <- objects$max_dist
      ncell_local <- (2 * max_dist + 1)^2 
      step_size  <- objects$step_size
      n_obs      <- length(objects$obs) + 1
      k_exp   <- ifelse(sim, exp(par[length(par)]), exp(par[length(par) - 1]))
      bg_rate <- ifelse(sim, 0, plogis(par[length(par)]))

      # Matrices of dimension (n_obs * ncell_local) rows x 9 columns
      env_i      <- objects$env_i
      nbhd_i     <- objects$nbhd_i
      to_dest    <- objects$to_dest
      dest       <- objects$dest

      # Build dispersal kernel
      kernel <- calculate_dispersal_kernel(max_dispersal_dist = step_size, 
                    kfun = function(x) dexp(x, k_exp))

      # Calculate environmental attraction
      attract0 <- env_function(env_i, par, nbhd = nbhd_i, sim = sim) 
      attract <- apply_kernel(attract0, kernel)

      # 3D array for propagating probabilities forward
      # At step 1, all probs initialized at 1 at the center cell of each nbhd
      current <- array(0, dim = c(ncell_local, n_obs, self$propagation_steps))
      center <- ncell_local / 2 + 0.5
      current[center, , 1] <- 1 

      # Propagation loop
      for (j in 1:(self$propagation_steps - 1)) {
        step_prob <- as.vector(current[, , j]) * attract[]
        dest[] <- step_prob[objects$to_dest_vec]
        p_step <- rowSums(dest, na.rm = TRUE)
        p_with_bg <- p_step + bg_rate - p_step * bg_rate %>%
          matrix(nrow = ncell_local, ncol = n_obs)
        col_totals <- colSums(p_with_bg)
        current[, , j + 1] <- p_with_bg / rep(col_totals, each = ncell_local)
      }
      return(current)
    },

    ## log_likelihood ----------------------------------------------------------
    log_likelihood = function(par, objects, sim, debug = FALSE) {

      obs        <- objects$obs
      outliers   <- objects$outliers
      n_obs      <- length(obs) + 1

      current <- self$build_kernel(par, objects, sim)

      # Calculate log likelihood 
      predictions <- matrix(0, nrow = self$propagation_steps, ncol = n_obs)
      for (i in 1:n_obs) {
        if (i %in% outliers || is.na(obs[i])) {
          predictions[, i] <- rep(NA, self$propagation_steps)
          next
        }
        predictions[, i] <- current[obs[i], i, ]
      }      
      log_likelihood <- rowSums(log(predictions), na.rm = TRUE) 
      out            <- -max(log_likelihood, na.rm = TRUE)
      print(paste("Current LL:", round(-out, 4)))
      print(paste("PP parameters:", paste(round(par, 4), collapse = ", ")))
      # if (is.infinite(out) || is.na(out)) out <- 0
      return(out)
    },

    ## fit ---------------------------------------------------------------------
    fit = function(trajectory, max_dist, step_size, rdf, par_start, sim = FALSE, 
                  outliers = integer(0), dt) {
                    
      objects <- self$prepare_objects(trajectory, max_dist, step_size, rdf, sim)
      objects$outliers <- outliers
      lbound <- if (sim) {
        c(rep(-Inf, length(par_start)))
      } else {
        c(rep(-Inf, length(par_start) - 2), -5, -30)
      }

      ubound <- if (sim) {
        c(rep(Inf, length(par_start)))
      } else {
        c(rep(Inf, length(par_start) - 2), 5, -2)
      }
      # Fit model
      tryCatch({
        par_out <- optim(par_start, self$log_likelihood, objects = objects,                         
                        sim = sim,  method = "L-BFGS-B",
                        lower = lbound,
                        upper = ubound,
                        control = list(
                          maxit = 2000,     # More iterations
                          factr = 1e5
                        ))
                        
        ll <- self$log_likelihood(par_out$par, objects, sim)
        
        return(list(
          par = par_out$par,
          ll = ll,
          convergence = par_out$convergence
          # objects = objects
        ))
        
      }, error = function(e) {
        message(paste("Path propagation fitting error:", e$message))
        return(NA)
      })
    },


    dispersal_from = function(init_point, par, step_size, n_steps, rdf = brdf,
                              threshold = 1e-6) {
      
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
          env_kernel  <- env_function(rdf[nbhd, 1:6], par, sim = FALSE)
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
    }
  )
)

## Theoretical simulations =====================================================

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
                                neighb = self$config$step_size,
                                rdf = rdf)

              # Plot path on landscape
              plot_pdf(nm = paste0("figs/sims/current", self$config$name, "_path_", i,
                  ".pdf"), x = 6, y = 6)
              terra::plot(land_i, xlim = c(100, 300), ylim = c(100, 300),
                main = paste0("Path ", i, " (", self$config$name, ")"))
              points(path_i$x, path_i$y, col = rgb(1, 1, 1, 0.3), pch = 19, cex = 0.5)
              lines(path_i$x, path_i$y, col = rgb(1, 1, 1, 0.3), lwd = 1.5)
              dev.off()
              
              results[[i]] <- list(path = path_i, rdf = rdf)
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
            
            # Build individual neighborhood
            nbhd_full <<- make_nbhd(rdf = rdf, i = seq_len(nrow(rdf)), 
                              sz = step_size)
            
            pp_model$propagation_steps <- self$config$propagation_steps()

            # Fit both models
            ss_result <- ss_model$fit(trajectory, max_dist, rdf = rdf, 
                        par_start, sim = TRUE)
            
            pp_result <- pp_model$fit(trajectory, max_dist, rdf = rdf,
              par_start, sim = TRUE, step_size = step_size)
            
            return(list(
              step_selection = ss_result,
              path_propagation = pp_result
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
              
              # Build individual neighborhood
              nbhd_full <<- make_nbhd(rdf = rdf, i = seq_len(nrow(rdf)), 
                                    sz = self$config$step_size)
              
              # Fit both models
              ss_result <- ss_model$fit(trajectory, max_dist, rdf = rdf, 
                         par_start = par_start, sim = TRUE)
              pp_result <- pp_model$fit(trajectory, max_dist, rdf = rdf, 
                         par_start = par_start, sim = TRUE, step_size = step_size)
              results[[i]] <- list(
                step_selection = ss_result,
                path_propagation = pp_result
              )
          }
        }
      
      return(results)
    }
  )
)

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
        step_size = NULL,
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
                        step_size = 1, obs_interval = 1, n_steps = 1000) {
            self$name <- name
            self$env_size <- env_size
            self$autocorr_range <- autocorr_range
            self$autocorr_strength <- autocorr_strength
            self$env_response <- env_response
            self$step_size <- step_size
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
        max_dist = function() self$step_size * (self$obs_interval + 1),
        propagation_steps = function() self$obs_interval + 2,

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
                obs_interval      = self$obs_interval
            )
            saveRDS(params, filepath)            
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
        paths <- sim$simulate_paths()                          # Generate paths
        fit_results <- sim$fit_models(paths, parallel, n_cores) # Fit models
        self$results[[config$name]] <- list(                    # Store results
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
      summary_list <- list()
      results_list <- list()
      
      for (name in names(self$results)) {
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
        ll_diff <- (ll_step - ll_pp)  # Log-likelihood difference
        ll_diff_per_obs <- ll_diff / (n_obs - 1)  # Per observation

        results_list[[name]] <- data.frame(
          individual = seq_along(ll_step),
          ll_step = ll_step,
          ll_pp = ll_pp
        )
        
        # Calculate AIC difference (positive means path propagation better)
        aic_diff <- 2 * (ll_step - ll_pp)
        
        summary_list[[name]] <- data.frame(
          config_name = name,
          autocorr_range = config$autocorr_range,
          barrier_density = config$b_density,
          obs_interval = config$obs_interval,
          mean_aic_diff = mean(aic_diff, na.rm = TRUE),
          median_aic_diff = median(aic_diff, na.rm = TRUE),
          ll_diff_per_obs = ll_diff_per_obs,
          prop_pp_better = mean(aic_diff > 0, na.rm = TRUE),
          n_successful = sum(!is.na(aic_diff))
        )
      }
      
      return(list(do.call(rbind, results_list), do.call(rbind, summary_list)))
    }
))

## Empirical data fitting ======================================================

# Jaguar class to handle data related to a single individual
#   $get_track() - get processed track data with movement metrics
#   $get_track_cells() - get raster cell indices for the track
#   $get_landscape() - get cropped landscape raster and dataframe for the track
jaguar <- R6Class("jaguar",
  public = list(
    id = NULL,

    initialize = function(id) {
      self$id <- as.numeric(id)
    },

    get_track = function() {
      dat <- jag_move[ID == self$id]
      dat$timestamp <- as.POSIXct(dat$timestamp, 
                              format = "%m/%d/%Y %H:%M")
      dat$year <- as.numeric(format(dat$timestamp, "%Y"))
      dat$year <- ifelse(dat$year > 23, dat$year + 1900, dat$year + 2000)

      dat$mon <- as.numeric(format(dat$timestamp, "%m"))
      dat$day <- as.numeric(format(dat$timestamp, "%d"))
      dat$hr <- format(dat$timestamp, "%H:%M")
      dat$hr <- as.numeric(gsub(":[0-9][0-9]", "", dat$hr))
      
      path <- jag_move[ID == self$id]
      path$t <- lubridate::mdy_hm(as.character(path$timestamp))
      path <- vect(path, geom = c("longitude", "latitude"), crs = wgs84)
      path <- project(path, epsg5880)
      path <- track(x = crds(path)[, 1], y = crds(path)[, 2], 
                    t = path$t, id = path$ID, crs = epsg5880)
      st <- steps(path)
      dat$sl <- c(NA, st$sl_)             # step lengths in m
      dat$ta <- c(NA, st$ta_)             # turn angles in radians
      dat$dir <- c(NA, st$direction_p)    # bearing in radians
      dat$dt <- c(NA, as.numeric(st$dt_)) # time interval in minutes
      dat$spd <- dat$sl / dat$dt
      return(dat[, c("timestamp", "longitude", "latitude", "ID", "year", "mon", 
                    "day", "hr", "sl", "ta", "dir", "dt", "spd")])
    },

    get_track_cells = function() {
      track <- self$get_track()
      return(cellFromXY(brazil_ras, 
                         track[, c("longitude", "latitude")]))
    },

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

    get_fragmentation = function() {
      # get maximum displacement radius from starting point
      track <- self$get_track()
      track_extent <- ext(terra::vect(track, 
                          geom = c("longitude", "latitude"), crs = wgs84))

      start_point <- c(track$longitude[1], track$latitude[1])
      start_cell <- cellFromXY(brazil_ras, matrix(start_point, nrow = 1))
      
    }
  ))

# Batch runner for empirical data fitting
# Needs a bit of refactoring for cyclomatic complexity
#   $run_all() - run fitting for all individuals in config
#   $process_individual() - process a single individual 
empirical_batch <- R6Class("empirical_batch",
  public = list(
    config = NULL,
    results = list(),

    initialize = function(config) {
      self$config <- config
    },

    run_all = function() {
      # Determine which individuals to process
      if (is.null(self$config$individuals)) {
        individuals_to_process <- jag_id$jag_id
      } else {
        individuals_to_process <- self$config$individuals
      }
      
      # Check for already completed results
      outfiles <- list.files("data/output")
      done <- gsub(".rds", "", outfiles) %>% 
        gsub("out_",     "", .) %>% 
        gsub("ll_null_", "", .) %>%
        gsub("ll_",      "", .) %>% 
        as.numeric()
      done <- done[!is.na(done)]
      i_todo <- setdiff(individuals_to_process, done)
      
      message(paste0("Processing ", length(i_todo), " individuals"))
      
      # Set up global neighborhood if needed
      if (!exists("nbhd_full")) {
        message("Generating global neighborhood matrix")
        nbhd_full <<- make_nbhd(i = seq_len(nrow(brdf)), sz = self$config$step_size)
      }
      
      # Set up parallel processing
      if (self$config$parallel) {
        registerDoParallel(self$config$n_cores)
        message(paste0("Using ", self$config$n_cores, " cores"))
      }
      
      # Create model instances
      ss_model <- step_selection_model$new()
      pp_model <- path_propagation_model$new()
      pp_model$propagation_steps <- self$config$propagation_steps 
      
      # Main processing loop
      if (self$config$parallel) {
        results <- foreach(i = i_todo, .packages = c("terra", "ctmm", "amt")) %dopar% {
          self$process_individual(i, ss_model, pp_model)
        }
      } else {
        results <- list()
        for (i in i_todo) {
          results[[length(results) + 1]] <- self$process_individual(i, ss_model, pp_model)
        }
      }
      
      # If there are already completed results, load and combine
      if (length(done) > 0) {
        existing_results <- list()
        for (i in done) {
          if (file.exists(paste0("data/output/out_", i, ".rds"))) {
            existing_results[[length(existing_results) + 1]] <- 
              readRDS(paste0("data/output/out_", i, ".rds"))
          } else {
            existing_results[[length(existing_results) + 1]] <- NA
          }
        }
        results <- c(existing_results, results)
      }
      self$results <- results

      # Clean up temporary output files 
      outfiles <- c(list.files("data/output", pattern = "^out_.*\\.rds$"),
                    list.files("data/output", pattern = "^NA_.*\\.rds$"))
      file.remove(file.path("data/output", outfiles))

      return(results)
    },
    
    process_individual = function(i, ss_model, pp_model) {
      message(paste0("Processing jaguar #", i))
      # Create jaguar instance and get data
      jag <- jaguar$new(i)
      track <- jag$get_track()
      track_cells <- jag$get_track_cells()
      n_obs <- length(track_cells)

      dt_scaled <- track$dt[2:length(track$dt)] / median(na.exclude(track$dt))
      dt_discrete <- round(dt_scaled)
      outliers <- which(dt_discrete != 1)

      if (length(outliers) > 0) dt_discrete <- dt_discrete[-outliers]
      
      if (self$config$holdout_set && nrow(track) > 100) {
        hold <- seq_len(ceiling(nrow(track) * self$config$holdout_frac))
        if (self$config$fit_model) {
          track_cells <- track_cells[hold]
        } else {
          track_cells <- track_cells[-hold]
          outliers <- outliers[outliers > max(hold)]
          if (length(outliers) > 0) outliers <- outliers - length(hold)
        }
      }
      
      # Calculate max distance for this individual
      pp_model$propagation_steps <- 8
      sl_emp <- na.exclude(track$sl[-outliers])
      max_dist <- ceiling(1.1 * max(sl_emp) / 1000)
            
      # Starting parameters
      par_start <- c(rep(0, self$config$npar - 2), # no environment preferences
                    log(1.0),                      # k_exp 
                    -15)                           # bg_rate ~ 0 

      # Fit models based on config
      if (self$config$model_type == 1) {
        # Step selection
        result <- ss_model$fit(track_cells, max_dist, rdf = brdf, 
                               par_start = par_start, outliers = outliers)
      } else if (self$config$model_type == 2) {
        # Path propagation  
        result <- pp_model$fit(track_cells, max_dist, self$config$step_size, 
                            rdf = brdf, par_start = par_start,
                            outliers = outliers)
      }
      
      # Save individual result
      if (!is.na(result[1])) {
        saveRDS(result, paste0("data/output/out_", i, ".rds"))
      } else {
        saveRDS(NA, paste0("data/output/NA_", i, ".rds"))
      }
      
      return(result)
    }
  ))

## Results and analysis ========================================================

# Builds on 'jaguar' to load fitted results and perform analyses
#    $load_results(): loads results from files
#    $compare_dispersal(): compares fitted dispersal kernels from both models
individual_analysis <- R6Class("individual_analysis",
  public = list(
    id = NULL,
    jaguar = NULL,
    track = NULL,
    track_cells = NULL,
    landscape = NULL,
    results = NULL,

    file_ss = NULL,
    file_pp = NULL,

    ss_disp = NULL,
    pp_disp = NULL,

    initialize = function(id = NULL, file_ss, file_pp) {
      self$id <- as.numeric(id)
      self$jaguar <- jaguar$new(id)
      self$track <- self$jaguar$get_track()
      self$track_cells <- self$jaguar$get_track_cells()
      self$landscape <- self$jaguar$get_landscape()
      self$file_ss <- file_ss
      self$file_pp <- file_pp
      self$results <- self$load_results(id = self$id)
    },

    load_results = function(id = NULL) {
      if (!is.null(id)) {
        self$id <- as.numeric(id)
      }
      all_results <- results_table(self$file_ss, self$file_pp)
      return(all_results[all_results$ID == self$id, ])
    },

    compare_dispersal = function(init_point = NULL, step = NULL, step_size, n_steps) {

      if (is.null(step)) step <- 1
      init_point <- self$track_cells[step]

      max_displacement <- step_size * n_steps

      fitted_pars_ss <- self$results[3:10]
      fitted_pars_pp <- self$results[14:21]

      ss_model <- path_propagation_model$new()
      pp_model <- path_propagation_model$new() # Going to unify classes dw

      self$ss_disp <- ss_model$dispersal_from(init_point = init_point, 
                        par = fitted_pars_ss, step_size = step_size, 
                        n_steps = n_steps)
      message("Step selection dispersal calculated")

      self$pp_disp <- pp_model$dispersal_from(init_point = init_point, 
                         par = fitted_pars_pp, step_size = step_size, 
                         n_steps = n_steps)
      message("Path propagation dispersal calculated")

      saveRDS(list(ss = self$ss_disp, pp = self$pp_disp), 
              paste0("data/output/dispersal_", self$id, "_", Sys.Date(), ".rds"))

      ########## PLOTTING STARTS HERE (probably move somewhere else) ###########

      plot_pdf(nm = paste0("figs/indiv_test_", Sys.Date(), ".pdf"), x = 8, y = 8)
      par(mfrow = c(1, 2))
      # Extract probability distributions
      p_ss <- self$ss_disp$probs
      p_pp <- self$pp_disp$probs

      # Create rasters using full computational area
      center <- (2 * max_displacement + 1)^2 / 2 + 1
      p_ss[center] <- NA  # Set center to NA for better color scaling
      p_pp[center] <- NA

      # Scale for visualization
      # p_ss_scaled <- sqrt(p_ss / max(p_ss, na.rm = TRUE))
      # p_pp_scaled <- sqrt(p_pp / max(p_pp, na.rm = TRUE))
      p_ss_scaled <- log(p_ss + 1) / log(max(p_ss, na.rm = TRUE) + 1)
      p_pp_scaled <- log(p_pp + 1) / log(max(p_pp, na.rm = TRUE) + 1)

      # Create rasters
      rast_ss <- rast(matrix(p_ss_scaled, nrow = 2 * max_displacement + 1, ncol = 2 * max_displacement + 1))
      rast_pp <- rast(matrix(p_pp_scaled, nrow = 2 * max_displacement + 1, ncol = 2 * max_displacement + 1))

      # Define plotting extent based on max_displacement``
      init_coords <- xyFromCell(brazil_ras, init_point)
      pixel_res <- res(brazil_ras)[1]
      half_pix <- max_displacement * pixel_res
      plot_extent <- ext(c(
        init_coords[1] - half_pix, init_coords[1] + half_pix,
        init_coords[2] - half_pix, init_coords[2] + half_pix
      ))

      crs(rast_ss) <- crs(brazil_ras)
      ext(rast_ss) <- plot_extent
      crs(rast_pp) <- crs(brazil_ras)
      ext(rast_pp) <- plot_extent

      # Plotting
      track_coords <- self$track[, c("longitude", "latitude")]
      terra::plot(rast_ss, main = "Step selection")
      points(track_coords, pch = 19, cex = 0.3, col = rgb(1, 0, 0, 0.3))
      terra::plot(rast_pp, main = "Path propagation")
      points(track_coords, pch = 19, cex = 0.3, col = rgb(1, 0, 0, 0.3))

      # # Difference plot
      # diff <- p_pp_scaled - p_ss_scaled
      # rast_diff <- rast(matrix(diff, nrow = 2 * max_displacement + 1, ncol = 2 * max_displacement + 1))
      # terra::plot(rast_diff, main = "Difference (PP - SS)")
      # # Forest cover context
      # forest_cover <- terra::crop(brazil_ras[[4]], plot_extent)
      # terra::plot(forest_cover, main = "Forest cover")
      dev.off()
      message("Dispersal comparison plot saved")                
    }
  ))