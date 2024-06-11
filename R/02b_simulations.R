# Simulation of animal movement and habitat selection, with parallel processing

source("R/00_functions.R")

## Switches ====================================================================

simname <- "s16all"

# Switches for reusing old data
gen_land   <- F
gen_path   <- F

# Switches for fitting models
fit_indiv  <- F
fit_all    <- T
debug_fit  <- F

# Number of cores to use for path generation and fitting
ncore_path <- 10
ncore_fit  <- 5

## Parameters ==================================================================

### Landscape generation parameters:
envsize <- 400    # size of landscape in cells
s1 <- 1           # strength of autocorrelation 
r1 <- 80          # range of autocorrelation in cells

### Model parameters:
# Order: par_move, par_env0, par_env1, par_env2
# par0 <- c(1, 3, -2, 0.3)        
par0 <- c(1, 2, -0.2, -0.2)

### Path generation parameters:
step_size    <- 1            # Max # pixels for each step
obs_interval <- 0             # Number of steps to skip between observations
n_step       <- 2000          # Number of steps to simulate
sim_n        <- 20           # Number of simulations 
n_obs        <- ceiling(n_step / (obs_interval + 1))

### Write parameters to file
params <- list(envsize, s1, r1, obs_interval, n_step, n_obs, sim_n, step_size, 
               par0[[1]], par0[[2]], par0[[3]], par0[[4]])
names(params) <- c("envsize", "s1", "r1", "obs_interval", "n_step", "n_obs", 
                   "sim_n", "step_size", "par_move", "par_env0", "par_env1", 
                   "par_env2")
saveRDS(params, "simulations/params.rds")
print(params)

if (any(is.na(par0))) par0 <- par0[!is.na(par0)]

# Value to start fitting from
par_start <- rep(1, length(par0))
# par_start <- par0

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
    # paths <- lapply(1:sim_n, function(i) {
    paths <- foreach(i = 1:sim_n, .packages = "terra") %dopar% {
        env01 <- unwrap(env02)
        message(paste0("Path #: ", i, " / ", sim_n))
        x0 <- ceiling(envsize / 2)
        y0 <- ceiling(envsize / 2)
        jag_path(x0, y0, n_step, par = par0, neighb = step_size)
        }
    # )
    saveRDS(paths, "simulations/paths.rds")
    message("Saved paths.")
    registerDoSEQ()
}

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


max_dist <- step_size * (obs_interval + 1)
ncell_local <- (2 * max_dist + 1) ^ 2
sim_steps   <- obs_interval * step_size + 2
# Number of steps to simulate, interval + first and last steps

## Fitting =====================================================================
if (fit_indiv || fit_all) {
    parallel_setup(ncore_fit)
    message("Fitting model parameters")
    done <- list.files("simulations", pattern = "par_out_")
    done <- gsub("par_out_", "", done) %>%
            gsub(".rds", "", .) %>%
            as.numeric()
    todo <- setdiff(1:sim_n, done)
    message(paste0("Fitting ", length(todo), " individuals"))
    env02 <- terra::wrap(env01) # foreach needs this

    if (fit_indiv) {
        # fit <- do.call(rbind, lapply(todo, function(i) {
        foreach(i = todo, .combine = rbind, .packages = "terra") %dopar% {
            message(paste0("Fitting individual #: ", i, " / ", length(todo)))

            env <- unwrap(env02)
            traject <- jag_traject_cells[[i]]
            prep_model_objects(traject, max_dist, env)
            
            message("Fitting parameters for model 1: path-dependent kernel")
            objects1 <- list(env, nbhd, max_dist, sim_steps, to_dest, obs)
            par_out1 <- optim(par_start, log_likelihood, objects = objects1)
            ll <- log_likelihood(par_out1$par, objects1)
            message(paste0("Saving log-likelihood for model 1: ", i))
            saveRDS(ll, file = paste0("simulations/ll_fit1", i, ".rds"))
            saveRDS(par_out1$par,
                    file = paste0("simulations/par_out_", i, ".rds"))
            message(paste0("COMPLETED path #: ", i, " / ", sim_n))
        }
        # ))    
    } else if (fit_all) {
        opt_fun <- function(par, env, traj) {
            env02 <- wrap(env)
            l <- foreach(i = 1:sim_n, .combine = c, .packages = "terra",
                         .export = ls(globalenv())) %dopar% {
                message("Fitting all individuals")
                env <- unwrap(env02)
                traject <- traj[[i]]
                prep_model_objects(traject, max_dist, env)

                objects1 <- list(env, nbhd, max_dist, sim_steps, to_dest, obs)
                ll <- log_likelihood(par, objects1)
                return(ll)
                }
            return(sum(l, na.rm = T))
        }
        par_out <- optim(par_start, opt_fun, env = env01, traj = jag_traject_cells)
        message("Completed fitting all individuals")
        saveRDS(par_out$par, file = "simulations/par_out_all.rds")
    }
}

## Cleanup =====================================================================

system(paste0("mkdir simulations/", simname))
system(paste0("mv simulations/*.rds simulations/", simname, "/."))
system(paste0("mv simulations/*.tif simulations/", simname, "/."))
system(paste0("cp simulations/", simname, "/env* simulations/."))

## Scratch =====================================================================

## Traditional SSF fitting code
    # message("Fitting parameters for model 2: traditional SSF") #------------ 
    # obs_points <- as.data.frame(jag_traject[[i]]) 
    # names(obs_points) <- c("x", "y")
    # tr <- amt::steps(amt::make_track(obs_points, x, y))
    # sl_emp <- as.vector(na.exclude(tr$sl_))
    # ta_emp <- as.vector(na.exclude(tr$ta_))
    # mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, 
    #                            max_dist = max_dist, scale = 1)
    # objects2 <- list(env, max_dist, mk, obs)
    # par_out2 <- optim(par_start, log_likelihood0, objects = objects2)
    # ll <- log_likelihood0(par_out2$par, objects2)
    # message(paste0("Saving log-likelihood for model 2: ", i))
    # saveRDS(ll, file = paste0("simulations/ll_fit2", i, ".rds"))
