# Simulation of animal movement and habitat selection, with parallel processing
source("R/00_functions.R")

## Switches ====================================================================

simname <- "s1"

# Switches for landscape generation, path generation, fitting, and plotting.
gen_land   <- F
gen_path   <- T
fit_indiv  <- T
plot_results <- F

# Roads once taken
fit_all    <- F
debug_fit  <- F
demo       <- F

# Number of cores to use for path generation and fitting
ncore_path <- 5
ncore_fit  <- 5

## Parameters ==================================================================

### Landscape generation parameters:
envsize <- 200    # size of landscape in cells
s1 <- 5           # strength of autocorrelation 
r1 <- 15          # range of autocorrelation in cells

### Model parameters:
# Order: par_env0, par_env1, par_env2, par_kexp, par_bgrate
par0 <- c(3, -2, 0.3)        
# par0 <- c(NA, 2, -0.2, -0.2)

### Path generation parameters:
step_size    <- 1             # Max # pixels for each step
obs_interval <- 0             # Number of steps to skip between observations
n_step       <- 1000          # Number of steps to simulate
sim_n        <- 20            # Number of simulations 
n_obs        <- ceiling(n_step / (obs_interval + 1))

sim_steps  <- obs_interval + 2 # Number of steps to simulate forward in PP model

### Write parameters to file
if (gen_land || gen_path || fit_indiv || fit_all) {
    params <- list(envsize, s1, r1, obs_interval, n_step, n_obs, sim_n, step_size, 
                par0[[1]], par0[[2]], par0[[3]])
    names(params) <- c("envsize", "s1", "r1", "obs_interval", "n_step", "n_obs", 
                    "sim_n", "step_size", "par_env0", "par_env1", "par_env2")
    saveRDS(params, "simulations/params.rds")
    print(params)
}

if (any(is.na(par0))) par0 <- par0[!is.na(par0)]

# Value to start fitting from
# par_start <- rep(0, length(par0))
par_start <- par0

## Landscape ===================================================================
if (!gen_land) {
    message("Reusing old landscape")
    env01 <- rast("simulations/env01.tif")
} else {
    message("Generating new landscape")
    env01 <- gen_landscape(size = envsize, s = s1, r = r1)
    terra::plot(env01)
    writeRaster(env01, "simulations/env01.tif", overwrite = TRUE)
}

## Paths =======================================================================
if (!gen_path) {
    message("Reusing old paths")
    paths <- readRDS("simulations/paths.rds")
} else {
    message("Simulating new paths")
    parallel_setup(ncore_path)
    env02 <- terra::wrap(env01) # foreach needs this
    paths <- foreach(i = 1:sim_n, .packages = c("metaRange", "terra")) %dopar% {
    # paths <- for (i in 1:sim_n) { # easier to debug
        env01 <- unwrap(env02)
        message(paste0("Path #: ", i, " / ", sim_n))
        x0 <- ceiling(envsize / 2)
        y0 <- ceiling(envsize / 2)
        return(jag_path(x0, y0, n_step, par = par0, neighb = step_size))
    }
    # )
    saveRDS(paths, "simulations/paths.rds")
    message("Saved paths.")
    registerDoSEQ()
}

## Fitting =====================================================================
if (fit_indiv || fit_all) {
    ### Prepare to fit 
    jag_traject <- lapply(paths, function(p) {
        out <- cbind(p$x, p$y)
        ind <- seq(1, nrow(out), obs_interval + 1)
        out <- out[ind, ]
        return(out)
    })
    jag_traject_cells <- lapply(jag_traject, function(tr) {
        out <- terra::cellFromRowCol(env01, tr[, 1], tr[, 2])
        return(out)
    })
    rdf <- raster_to_df(env01)
    # Make global neighborhood from raster
    message("Building global neighborhood")
    nbhd0 <- make_nbhd(i = seq_len(ncell(env01)), sz = step_size, r = env01, rdf = rdf) 

    max_dist <- step_size * (obs_interval + 1)
    ncell_local <- (2 * max_dist + 1) ^ 2
    # sim_steps   <- obs_interval * step_size + 2
    # Number of steps to simulate, interval + first and last steps

    if (!gen_path) parallel_setup(ncore_fit)

    message("Fitting model parameters")
    done <- list.files("simulations", pattern = "par_out_")
    done <- gsub("par_out_", "", done) %>%
            gsub(".rds", "", .) %>%
            as.numeric()
    todo <- setdiff(1:sim_n, done)
    message(paste0("Fitting ", length(todo), " individuals"))

    # Prepare objects for fitting
    objects_all <- lapply(seq_len(sim_n), function(n) {
        paste0("Preparing objects for path: ", n) %>% message()
        out           <- prep_model_objects(jag_traject_cells[[n]], max_dist, 
                                            rdf = rdf, sim = TRUE)
        out$sim_steps <- sim_steps
        return(out)
    })
    
    if (fit_indiv) {
        # Fit individuals one at a time ----------------------------------------
        # fit <- do.call(rbind, lapply(todo, function(i) { # easier to debug
        foreach(i = todo, .combine = rbind, .packages = c("terra", "metaRange")) %dopar% {
            message(paste0("Fitting individual #: ", i, " / ", length(todo)))            
            message("Fitting parameters for model 1: step-selection")
            objects1 <- objects_all[[i]]
            
            par_out1 <- optim(par_start, log_likelihood0, objects = objects1)
            ll1 <- log_likelihood0(par_out1$par, objects1)
            message("Fitting parameters for model 2: path-propagation")
            par_out2 <- optim(par_start, log_likelihood, objects = objects1)
            ll2 <- log_likelihood(par_out2$par, objects1)

            message(paste0("Saving log-likelihoods: ", i))
            saveRDS(ll1, file = paste0("simulations/ll_fit1", i, ".rds"))
            saveRDS(par_out1$par,
                    file = paste0("simulations/par_out1_", i, ".rds"))
            saveRDS(ll2, file = paste0("simulations/ll_fit2", i, ".rds"))
            saveRDS(par_out2$par,
                    file = paste0("simulations/par_out2_", i, ".rds"))
            message(paste0("COMPLETED path #: ", i, " / ", sim_n))
        }
        # ))    # easier to debug
    }
    
    if (fit_all) {
        # Fit all individuals at the same time ---------------------------------
        message("Fitting all individuals")
        # Exports need to be passed to foreach inside optim
        optim_function <- function(par) {
            l <- foreach(i = 1:sim_n, .combine = c, .packages = "terra",
                         .export = ls(globalenv())) %dopar% {
                    objects1 <- objects_all[[i]]
                    ll <- log_likelihood(par, objects1)
                    message(paste0("Log likelihood calculated, #: ", i))
                    return(ll)
                  }
            return(sum(l, na.rm = T))
        }

        par_out <- optim(par_start, optim_function)
        message("Completed fitting all individuals")
        saveRDS(par_out$par, file = "simulations/par_out_all.rds")
    }
}

## Cleanup =====================================================================

if (fit_indiv || fit_all) {
    system(paste0("mkdir simulations/", simname))
    system(paste0("mv simulations/*.rds simulations/", simname, "/."))
    system(paste0("mv simulations/*.tif simulations/", simname, "/."))
    system(paste0("cp simulations/", simname, "/env* simulations/."))
    system(paste0("cp simulations/", simname, "/paths.rds simulations/."))
    system(paste0("cp simulations/", simname, "/params.rds simulations/."))
}

## Demo ========================================================================

if (plot_results) {
    source("R/02c_simulation_analysis.R")
}
