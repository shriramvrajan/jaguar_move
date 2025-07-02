# Attempt to refactor into object-oriented idiom
library(R6)

## Landscape object for environmental raster
landscape <- R6Class("landscape", 
    public = list(
        raster = NULL,
        config = NULL,

        initialize = function(config) {
            self$config <- config
            message(paste("Generating landscape for", config$name))
            self$raster <- gen_landscape(
                size = config$env_size,
                s = config$autocorr_strength,
                r = config$autocorr_range
            )
        },

        save_raster = function(filepath) {
            writeRaster(self$raster, filepath, overwrite = TRUE)
        }
    )
)

# Step Selection Model Class
step_selection_model <- R6Class("step_selection_model",
  public = list(
    
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE, env_raster = NULL) {
      message("Preparing step selection objects")
      
      # Get environment data
      if (sim) {
        env_0 <- scale(rdf$sim1[])
        env_ras <- env_raster
      } else {
        env_0 <- scale(rdf[, 1])  # First column of brdf
        env_ras <- brazil_ras
      }
      if (any(is.na(env_0))) env_0[which(is.na(env_0))] <- 0
      
      # Extended neighborhoods for step selection
      nbhd_0 <- make_nbhd(i = trajectory, sz = max_dist, r = env_ras, rdf = rdf)
      
      # Build observations for step selection
      obs_0 <- sapply(seq_along(trajectory), function(i) {
        if (i == length(trajectory)) return(NULL)
        step <- trajectory[i + 1]
        return(which(nbhd_0[i, ] == step))
      }) %>% unlist()
      
      # Return step selection specific objects
      list(
        env_0 = env_0,
        nbhd_0 = nbhd_0, 
        obs_0 = obs_0,
        max_dist = max_dist,
        outliers = integer(0)  # Will be set externally if needed
      )
    },
    
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, env_raster = NULL, outliers = integer(0)) {
      # Prepare objects
      objects <- self$prepare_objects(trajectory, max_dist, rdf, sim, env_raster)
      objects$outliers <- outliers
      
      # Fit model
      tryCatch({
        par_out <- optim(par_start, log_likelihood0, objects = objects)
        ll <- log_likelihood0(par_out$par, objects)
        
        return(list(
          par = par_out$par,
          ll = ll,
          convergence = par_out$convergence,
          objects = objects
        ))
      }, error = function(e) {
        message(paste("Step selection fitting error:", e$message))
        return(NA)
      })
    }
  )
)

# Path Propagation Model Class  
path_propagation_model <- R6Class("path_propagation_model",
  public = list(
    
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE, env_raster = NULL, sim_steps = NULL) {
      message("Preparing path propagation objects")
      
      # Get environment data
      if (sim) {
        env_0 <- scale(rdf$sim1[])
        env_ras <- env_raster
      } else {
        env_0 <- scale(rdf[, 1])
        env_ras <- brazil_ras
      }
      if (any(is.na(env_0))) env_0[which(is.na(env_0))] <- 0
      
      # Extended neighborhoods
      nbhd_index <- make_nbhd(i = trajectory, sz = max_dist, r = env_ras, rdf = rdf)
      
      # Build observations for path propagation
      obs_i <- sapply(2:nrow(nbhd_index), function(r) {
        local <- nbhd_index[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })
      
      # Build complex nested neighborhood structure
      nbhd_list <- lapply(seq_len(nrow(nbhd_index)), function(i) {                
        row_inds <- seq_len(ncol(nbhd_index)) + (i - 1) * ncol(nbhd_index)
        names(row_inds) <- nbhd_index[i, ]
        out <- matrix(row_inds[as.character(nbhd0[nbhd_index[i, ], ])], 
                      nrow = length(row_inds), ncol = ncol(nbhd0))
        return(out)
      })
      nbhd_i <- do.call(rbind, nbhd_list)
      
      # Flatten nbhd_index for environment indexing
      nbhd_index <- as.vector(t(nbhd_index))
      env_i <- env_0[nbhd_index]
      if (any(is.na(env_i))) env_i[which(is.na(env_i))] <- 0
      
      # Build connectivity matrix
      to_dest <- tapply(seq_len(length(nbhd_i)), nbhd_i, function(x) {  
        if (length(x) <= ncol(nbhd_i)) {
          out <- c(x, rep(NA, ncol(nbhd_i) - length(x)))
        } else {
          out <- x[seq_len(ncol(nbhd_i))]
        }
        return(out)
      })
      to_dest <- t(matrix(unlist(to_dest), nrow = ncol(nbhd_i), ncol = nrow(nbhd_i)))
      dest <- matrix(0, nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
      
      # Return path propagation specific objects
      result <- list(
        env_i = env_i,
        nbhd_i = nbhd_i,
        to_dest = to_dest,
        dest = dest,
        obs_i = obs_i,
        max_dist = max_dist,
        outliers = integer(0)  # Will be set externally if needed
      )
      
      if (!is.null(sim_steps)) {
        result$sim_steps <- sim_steps
      }
      
      if (sim) {
        result$mu_env <- attributes(env_i)[[2]]
        result$sd_env <- attributes(env_i)[[3]]
      }
      
      return(result)
    },
    
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, env_raster = NULL, sim_steps = NULL, outliers = integer(0)) {
      # Prepare objects
      objects <- self$prepare_objects(trajectory, max_dist, rdf, sim, env_raster, sim_steps)
      objects$outliers <- outliers
      
      # Fit model
      tryCatch({
        par_out <- optim(par_start, log_likelihood, objects = objects)
        ll <- log_likelihood(par_out$par, objects)
        
        return(list(
          par = par_out$par,
          ll = ll,
          convergence = par_out$convergence,
          objects = objects
        ))
      }, error = function(e) {
        message(paste("Path propagation fitting error:", e$message))
        return(NA)
      })
    }
  )
)

## Movement simulator for path generation and model fitting 
movement_simulator <- R6Class("movement_simulator",
    public = list(
        config = NULL,

        initialize = function(config) {
            self$config <- config
        },

        simulate_paths = function(landscape, parallel = TRUE, n_cores = 5) {
            message(paste("Simulating", self$config$n_individuals, 
                            "paths for", self$config$name))
            
            if (parallel) {
                registerDoParallel(n_cores)
                env_wrapped <- terra::wrap(landscape$raster)
                
                paths <- foreach(i = 1:self$config$n_individuals, 
                                .packages = c("metaRange", "terra"),
                                .export = c("jag_path", "make_nbhd", "raster_to_df")) %dopar% {
                    env_unwrapped <- unwrap(env_wrapped)
                    message(paste0("Path #: ", i, " / ", 
                                   self$config$n_individuals))
                    x0 <- ceiling(self$config$env_size / 2)
                    y0 <- ceiling(self$config$env_size / 2)
                    return(jag_path(x0, y0, self$config$n_steps, 
                                    par = self$config$env_response, 
                                    neighb = self$config$step_size,
                                    env_raster = env_unwrapped))
                }
                registerDoSEQ()
            } else {
                paths <- list()
                for (i in 1:self$config$n_individuals) {
                message(paste0("Path #: ", i, " / ", self$config$n_individuals))
                x0 <- ceiling(self$config$env_size / 2)
                y0 <- ceiling(self$config$env_size / 2)
                paths[[i]] <- jag_path(x0, y0, self$config$n_steps, 
                                        par = self$config$env_response, 
                                        neighb = self$config$step_size,
                                        env_raster = landscape$raster)
                }
            }
            
            return(paths)
    },
    
    fit_models = function(paths, landscape, parallel = TRUE, n_cores = 5) {
        message(paste("Fitting models for", self$config$name))

        # Prepare trajectories
        jag_traject <- lapply(paths, function(p) {
            out <- cbind(p$x, p$y)
            ind <- seq(1, nrow(out), self$config$obs_interval + 1)
            out <- out[ind, ]
            return(out)
        })
        
        jag_traject_cells <- lapply(jag_traject, function(tr) {
            out <- terra::cellFromRowCol(landscape$raster, tr[, 1], tr[, 2])
            return(out)
        })
        
        rdf <- raster_to_df(landscape$raster)
        
        # Build global neighborhood
        message("Building global neighborhood")
        nbhd0 <<- make_nbhd(i = seq_len(ncell(landscape$raster)), 
                            sz = self$config$step_size, 
                            r = landscape$raster, rdf = rdf) 
        
        max_dist <- self$config$max_dist()
        sim_steps <- self$config$sim_steps()
        
        # Fit models
        par_start <- c(-1.5, 1.5, -0.2) # true value

        # Create model instances
        ss_model <- step_selection_model$new()
        pp_model <- path_propagation_model$new()
        
        results <- list()
        for (i in 1:self$config$n_individuals) {
            message(paste0("Fitting individual #: ", i))
            
            trajectory <- jag_traject_cells[[i]]
            
            # Fit both models
            ss_result <- ss_model$fit(trajectory, max_dist, rdf, par_start, 
                                    sim = TRUE, env_raster = landscape$raster)
            pp_result <- pp_model$fit(trajectory, max_dist, rdf, par_start,
                                    sim = TRUE, env_raster = landscape$raster, 
                                    sim_steps = sim_steps)
            
            results[[i]] <- list(
            step_selection = ss_result,
            path_propagation = pp_result
            )
        }
      
      return(results)
    }
  )
)

## Configuration parameters
simulation_config <- R6Class("simulation_config",
    public <- list(
        # Landscape parameters
        env_size = 200,
        autocorr_strength = 5,
        autocorr_range = 50,

        # Model parameters
        env_response = c(-1.5, 1.5, -0.2),
        step_size = 1,
        obs_interval = 5,

        # Simulation parameters
        n_steps = 1000,
        n_individuals = 10,

        # Output parameters
        name = "default", 

        initialize = function(name = "default", autocorr_range = 50, 
                            n_individuals = 10, env_size = 200, 
                            autocorr_strength = 5, 
                            env_response = c(-1.5, 1.5, -0.2),
                            step_size = 1, obs_interval = 5, n_steps = 1000) {
            self$name <- name
            self$autocorr_range <- autocorr_range
            self$n_individuals <- n_individuals
            self$env_size <- env_size
            self$autocorr_strength <- autocorr_strength
            self$env_response <- env_response
            self$step_size <- step_size
            self$obs_interval <- obs_interval
            self$n_steps <- n_steps
        },

        # Derived properties
        max_dist = function() max(2, self$step_size * (self$obs_interval + 1)),
        sim_steps = function() self$obs_interval + 2,

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

## Batch runner to coordinate multiple simulations
simulation_batch <- R6Class("simulation_batch",
  public = list(
    configs = list(),
    results = list(),
    
    # Create batch for fragmentation study
    autocorr_range_study = function(r1_values = c(10, 25, 50, 80), 
                                         n_individuals = 10) {
      self$configs <- lapply(r1_values, function(r1) {
        simulation_config$new(
          name = paste0("r1_", r1),
          autocorr_range = r1,
          n_individuals = n_individuals
        )
      })
      names(self$configs) <- sapply(self$configs, function(x) x$name)
      return(self)
    },
    
    # Run all simulations in the batch
    run_all = function(parallel = TRUE, n_cores = 5) {
      self$results <- list()
      
      for (config in self$configs) {
        message(paste("=" %>% rep(50) %>% paste(collapse = "")))
        message(paste("RUNNING SIMULATION:", config$name))
        message(paste("=" %>% rep(50) %>% paste(collapse = "")))
        
        # Create landscape
        land <- landscape$new(config)
        
        # Create simulator
        sim <- movement_simulator$new(config)
        
        # Generate paths
        paths <- sim$simulate_paths(land, parallel, n_cores)
        
        # Fit models
        fit_results <- sim$fit_models(paths, land, parallel, n_cores)
        
        # Store results
        self$results[[config$name]] <- list(
          config = config,
          landscape = land,
          paths = paths,
          fits = fit_results
        )
        
        message(paste("COMPLETED:", config$name))
      }
      
      return(self$results)
    },
    
    # Extract performance comparison for plotting
    get_performance_summary = function() {
      summary_list <- list()
      
      for (name in names(self$results)) {
        fits <- self$results[[name]]$fits
        config <- self$results[[name]]$config
        
        # Extract log-likelihoods for each individual
        ll_step <- sapply(fits, function(x) if(is.na(x[1])) NA else x$step_selection$ll)
        ll_pp <- sapply(fits, function(x) if(is.na(x[1])) NA else x$path_propagation$ll)
        
        # Calculate AIC difference (positive means path propagation better)
        aic_diff <- 2 * (ll_step - ll_pp)
        
        summary_list[[name]] <- data.frame(
          config_name = name,
          autocorr_range = config$autocorr_range,
          mean_aic_diff = mean(aic_diff, na.rm = TRUE),
          median_aic_diff = median(aic_diff, na.rm = TRUE),
          prop_pp_better = mean(aic_diff > 0, na.rm = TRUE),
          n_successful = sum(!is.na(aic_diff))
        )
      }
      
      return(do.call(rbind, summary_list))
    },
    
    # Quick plot of results
    plot_fragmentation_effect = function() {
      summary_df <- self$get_performance_summary()
      
      plotpdf(nm = paste0("figs/fragmentation_study_", Sys.Date(), ".pdf"), 
             x = 8, y = 6)
      par(mfrow = c(1, 2))
      
      # Plot 1: AIC difference vs fragmentation
      plot(summary_df$autocorr_range, summary_df$mean_aic_diff,
           xlab = "Autocorrelation Range (r1)", 
           ylab = "Mean AIC Difference (PP - SS)",
           main = "Path Propagation Advantage vs Fragmentation",
           pch = 19, cex = 1.5)
      abline(h = 0, lty = 2)
      text(summary_df$autocorr_range, summary_df$mean_aic_diff,
           labels = summary_df$config_name, pos = 3)
      
      # Plot 2: Proportion where PP is better
      plot(summary_df$autocorr_range, summary_df$prop_pp_better,
           xlab = "Autocorrelation Range (r1)",
           ylab = "Proportion where PP outperforms SS", 
           main = "Frequency of PP Advantage",
           pch = 19, cex = 1.5, ylim = c(0, 1))
      abline(h = 0.5, lty = 2)
      text(summary_df$autocorr_range, summary_df$prop_pp_better,
           labels = summary_df$config_name, pos = 3)
      
      dev.off()
      
      return(summary_df)
    }
  )
)
