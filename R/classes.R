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
    
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE,
                               env_raster = NULL) {
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
      obs <- sapply(2:nrow(nbhd_0), function(r) {
        local <- nbhd_0[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })
      
      # Return step selection specific objects
      list(
        env_0 = env_0,
        nbhd_0 = nbhd_0, 
        obs = obs,
        max_dist = max_dist,
        outliers = integer(0)  # Will be set externally if needed
      )
    },
    
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, 
                   env_raster = NULL, outliers = integer(0)) {

      # Prepare objects
      objects <- self$prepare_objects(trajectory, max_dist, rdf, 
                                      sim, env_raster)
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
    
    prepare_objects = function(trajectory, max_dist, rdf, sim = FALSE,
                               env_raster = NULL, sim_steps = NULL) {
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
      nbhd_0 <- make_nbhd(i = trajectory, sz = max_dist, r = env_ras, rdf = rdf)
      
      # Build observations for path propagation
      obs <- sapply(2:nrow(nbhd_0), function(r) {
        local <- nbhd_0[(r - 1), ]
        nextcell <- which(local == trajectory[r])
        if (length(nextcell) == 0) return(NA) else return(nextcell)
      })
      
      # Build complex nested neighborhood structure
      nbhd_list <- lapply(seq_len(nrow(nbhd_0)), function(i) {                
        row_inds <- seq_len(ncol(nbhd_0)) + (i - 1) * ncol(nbhd_0)
        names(row_inds) <- nbhd_0[i, ]
        out <- matrix(row_inds[as.character(nbhd_full[nbhd_0[i, ], ])], 
                      nrow = length(row_inds), ncol = ncol(nbhd_full))
        return(out)
      })
      nbhd_i <- do.call(rbind, nbhd_list)
      
      # Flatten nbhd_0 for environment indexing
      env_i <- env_0[as.vector(t(nbhd_0))]
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
      to_dest <- t(matrix(unlist(to_dest), 
                   nrow = ncol(nbhd_i), ncol = nrow(nbhd_i)))
      dest <- matrix(0, nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
      
      # Return path propagation specific objects
      result <- list(
        env_i = env_i,
        nbhd_i = nbhd_i,
        to_dest = to_dest,
        dest = dest,
        obs = obs,
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
    
    fit = function(trajectory, max_dist, rdf, par_start, sim = FALSE, 
                   env_raster = NULL, sim_steps = NULL, outliers = integer(0)) {
      # Prepare objects
      objects <- self$prepare_objects(trajectory, max_dist, rdf, sim, 
                                      env_raster, sim_steps)
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

        simulate_paths = function() {  # terra not working well with parallel
            message(paste("Simulating", self$config$n_individuals, 
                            "paths for", self$config$name))
            
            # Extract config values
            n_individuals <- self$config$n_individuals
            env_size <- self$config$env_size
            n_steps <- self$config$n_steps
            env_response <- self$config$env_response
            step_size <- self$config$step_size
            autocorr_strength <- self$config$autocorr_strength
            autocorr_range <- self$config$autocorr_range
        
            results <- list()
            for (i in 1:self$config$n_individuals) {
              set.seed(i + 5001)
              land_i <- gen_landscape(
                  size = self$config$env_size,
                  s = self$config$autocorr_strength,
                  r = self$config$autocorr_range
              )
              
              message(paste0("Path #: ", i, " / ", self$config$n_individuals))
              x0 <- ceiling(self$config$env_size / 2)
              y0 <- ceiling(self$config$env_size / 2)
              
              path_i <- jag_path(x0, y0, self$config$n_steps, 
                                par = self$config$env_response, 
                                neighb = self$config$step_size,
                                env_raster = land_i)
              
              results[[i]] <- list(path = path_i)
            }            
            return(results)
    },
    
    fit_models = function(paths, parallel = TRUE, n_cores = NULL) {
        message(paste("Fitting models for", self$config$name))
        
        # Extract config information
        n_cores <- n_cores %||% self$config$n_cores
        max_dist <- self$config$max_dist()
        sim_steps <- self$config$sim_steps()
        par_start <- c(-1.5, 1.5, -0.2)

        obs_interval <- self$config$obs_interval
        step_size <- self$config$step_size
        n_individuals <- self$config$n_individuals
        env_size <- self$config$env_size
        autocorr_strength <- self$config$autocorr_strength
        autocorr_range <- self$config$autocorr_range

        if (parallel) {
          registerDoParallel(n_cores) 
          # Get ALL functions from global environment
          all_functions <- ls(globalenv())[sapply(ls(globalenv()), function(x) is.function(get(x)))]
          
          results <- foreach(i = seq_len(n_individuals),
              .packages = c("terra", "gstat"),
              .export = c(  # Model classes
                            "step_selection_model", "path_propagation_model",
                            # All functions
                            all_functions,
                            # Data
                            "paths",
                            # Parameters
                            "max_dist", "sim_steps", "par_start", "obs_interval", 
                            "step_size", "env_size", "autocorr_strength", 
                            "autocorr_range")) %dopar% {                        
                    message(paste0("Fitting individual #: ", i))
      
                    path_i <- paths[[i]]$path
                    # Regenerate landscape deterministically
                    set.seed(i + 5001)
                    land_i <- gen_landscape(
                        size = env_size,
                        s = autocorr_strength,
                        r = autocorr_range
                    )                  

                    rdf <- raster_to_df(land_i)
                    trajectory_df <- cbind(path_i$x, path_i$y)
                    ind <- seq(1, nrow(trajectory_df), obs_interval + 1)
                    trajectory_df <- trajectory_df[ind, ]
                    trajectory <- terra::cellFromRowCol(land_i, 
                                    trajectory_df[, 1], trajectory_df[, 2])
                    
                    # Build individual neighborhood
                    nbhd_full <<- make_nbhd(i = seq_len(ncell(land_i)), 
                                      sz = step_size, r = land_i, rdf = rdf)
                    
                    ss_model <- step_selection_model$new()
                    pp_model <- path_propagation_model$new()
                    
                    # Fit both models
                    ss_result <- ss_model$fit(trajectory, max_dist, rdf, 
                                par_start, sim = TRUE, env_raster = land_i)
                    pp_result <- pp_model$fit(trajectory, max_dist, rdf,
                                par_start, sim = TRUE, env_raster = land_i, 
                                sim_steps = sim_steps)
                    
                    return(list(
                      step_selection = ss_result,
                      path_propagation = pp_result,
                      landscape = land_i
                    ))
                }
        } else {
          # Create model instances
          ss_model <- step_selection_model$new()
          pp_model <- path_propagation_model$new()
          
          results <- list()
          for (i in 1:self$config$n_individuals) {
              message(paste0("Fitting individual #: ", i))
              
              path_i <- paths[[i]]$path
              # Regenerate landscape deterministically
              set.seed(i + 5001)
              land_i <- gen_landscape(
                  size = env_size,
                  s = autocorr_strength,
                  r = autocorr_range
              )   
              
              rdf <- raster_to_df(land_i)
              # Prepare trajectory
              trajectory_df <- cbind(path_i$x, path_i$y)
              ind <- seq(1, nrow(trajectory_df), self$config$obs_interval + 1)
              trajectory_df <- trajectory_df[ind, ]
              trajectory <- terra::cellFromRowCol(land_i, 
                                    trajectory_df[, 1], trajectory_df[, 2])
              
              # Build individual neighborhood
              nbhd_full <<- make_nbhd(i = seq_len(ncell(land_i)), 
                                    sz = self$config$step_size, 
                                    r = land_i, rdf = rdf)
              # Fit both models
              ss_result <- ss_model$fit(trajectory, max_dist, rdf, par_start, 
                                      sim = TRUE, env_raster = land_i)
              pp_result <- pp_model$fit(trajectory, max_dist, rdf, par_start,
                                      sim = TRUE, env_raster = land_i, 
                                      sim_steps = sim_steps)
              results[[i]] <- list(
                step_selection = ss_result,
                path_propagation = pp_result,
                landscape = land_i
              )
          }
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
        obs_interval = 1,

        # Simulation parameters
        n_steps = 1000,
        n_individuals = 10,
        n_cores = 15,

        # Output parameters
        name = "default", 

        initialize = function(name = "default", autocorr_range = 50, 
                            n_individuals = 10, env_size = 200, 
                            autocorr_strength = 5, n_cores = 15,
                            env_response = c(-1.5, 1.5, -0.2),
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
        },

        # Derived properties
        # max_dist = function() max(2, self$step_size * (self$obs_interval + 1)),
        max_dist = function() self$step_size * (self$obs_interval + 1),
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
    autocorr_range_study = function(r1_values = c(1, 2, 4, 6, 8, 10, 15, 20, 30, 40, 50, 60, 70), 
                                         n_individuals = 30) {
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
    run_all = function(parallel = TRUE, n_cores = NULL) {
      self$results <- list()
      
      for (config in self$configs) {
        message(paste("=" %>% rep(50) %>% paste(collapse = "")))
        message(paste("RUNNING SIMULATION:", config$name))
        message(paste("=" %>% rep(50) %>% paste(collapse = "")))
        
        # Create simulator
        sim <- movement_simulator$new(config)
        n_cores <- n_cores %||% config$n_cores

        # Generate paths
        paths <- sim$simulate_paths()
        
        # Fit models
        fit_results <- sim$fit_models(paths, parallel, n_cores)
        
        # Store results
        self$results[[config$name]] <- list(
            config = config,
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
        ll_step <- sapply(fits, 
                        function(x) if (is.na(x[1])) NA else x$step_selection$ll)
        ll_pp <- sapply(fits, 
                      function(x) if (is.na(x[1])) NA else x$path_propagation$ll)
        
        # Calculate AIC difference (positive means path propagation better)
        aic_diff <- 2 * (ll_step - ll_pp)
        
        summary_list[[name]] <- data.frame(
          config_name = name,
          autocorr_range = config$autocorr_range,
          mean_aic_diff = mean(aic_diff, na.rm = TRUE),
          median_aic_diff = median(aic_diff, na.rm = TRUE),
          aic_diff = list(aic_diff),
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
             x = 4, y = 4)
      par(mfrow = c(1, 1))
      
      # Plot 1: AIC difference vs fragmentation
      plot(summary_df$autocorr_range, summary_df$median_aic_diff,
           xlab = "Autocorrelation Range (r1)", 
           ylab = "Median AIC Difference (PP - SS)",
           main = "Path Propagation Advantage vs Fragmentation",
           pch = 19, cex = 1.5)
      abline(h = 0, lty = 2)
      text(summary_df$autocorr_range, summary_df$median_aic_diff,
           labels = summary_df$config_name, pos = 3)
      
      dev.off()
      
      return(summary_df)
    }
  )
)
