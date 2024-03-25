# Simulation of animal movement and habitat selection, with parallel processing

source("R/00_functions.R")

## Switches ====================================================================

parallel_setup(20) # How many cores?
gen_land   <- F
gen_path   <- T
fit_indivs <- T

# take sim name and make a folder and save all the outputs there

## Parameters ==================================================================

### Landscape generation
envsize <- 200  # size of landscape in cells
s1 <- 0.01        # strength of autocorrelation 
r1 <- 0.01        # range of autocorrelation in cells

### Model parameters: env1 attraction scalar, move probability exponent
# par0   <- c(3, -2)    
# MOVE probability exponent, env1 attraction parameters 1 to 3
par0 <- c(2, 2, -0.2, -0.2)        

sim_interval <- 5             # GPS observations taken every n steps, for sim
n_step       <- 4000          # Number of steps to simulate
sim_n        <- 20           # Number of simulations 
step_size    <- 1             # Max # pixels for each step
n_obs        <- floor(n_step / sim_interval)

params <- list(envsize, s1, r1, sim_interval, n_step, n_obs, sim_n, step_size, 
               par0[[1]], par0[[2]], par0[[3]], par0[[4]])
names(params) <- c("envsize", "s1", "r1", "sim_interval", "n_step", "n_obs", 
                   "sim_n", "step_size", "par_move", "par_env0", "par_env1", 
                   "par_env2")
saveRDS(params, "data/output/simulations/params.rds")
print(params)

if (any(is.na(par0))) {
    par0 <- par0[-which(is.na(par0))]  # Remove NA values
}

## Landscape generation ========================================================

if (!gen_land) {
    message("Reusing old landscape")
    env01 <- list(rast("data/output/simulations/env01.tif"),
                  readRDS("data/output/simulations/env01.rds"))
} else {
    message("Generating new landscape")
    env01 <- gen_landscape(size = envsize, s = s1, r = r1)
    writeRaster(env01[[1]], "data/output/simulations/env01.tif", overwrite = TRUE)
    saveRDS(env01[[2]], "data/output/simulations/env01.rds")
    terra::plot(env01[[1]])
}

## Simulation ==================================================================

if (!gen_path) {
    message("Reusing old paths")
    paths <- readRDS("data/output/simulations/paths.rds")
} else {
    message("Simulating new paths")
    env02 <- terra::wrap(env01[[1]])
    # paths <- lapply(1:sim_n, function(i) {
    paths <- foreach(i = 1:sim_n, .packages = "terra") %dopar% {
        env01 <- list(unwrap(env02), env01[[2]])
        message(paste0("Path #: ", i, " / ", sim_n))
        x0 <- ceiling(envsize / 2)
        y0 <- ceiling(envsize / 2)
        jag_path(x0, y0, n_step, par = par0, neighb = step_size)
        }
    # )
    saveRDS(paths, "data/output/simulations/paths.rds")
    message("Saved paths.")
}

## HOW TO USE RASTER WITH FOREACH?

par(mfrow = c(4, 5))
for (i in 1:sim_n) {
    plot_path(paths[[i]])
}

## Fitting =====================================================================

if (fit_indivs) {
    sim_steps    <- sim_interval * 2  # Number of steps to simulate

    envdf <- env01[[2]]
    env_index <- seq_len(nrow(envdf))
    envdf <- cbind(envdf, env_index)

    jag_traject <- lapply(paths, function(p) {
        out <- cbind(p$x, p$y)
        ind <- seq(1, nrow(out), sim_interval)
        out <- out[ind, ]
        return(out)
    })

    jag_traject_cells <- lapply(jag_traject, function(tr) {
        raster::cellFromXY(env01[[1]], tr[, 1:2])
    })

    dist <- lapply(jag_traject, function(tr) {
        out <- c(0, sqrt(diff(tr[, 1])^2 + diff(tr[, 2])^2))
        return(out)
    })
    max_dist <- ceiling(max(unlist(dist)) * 1.5)
    step_range <- (2 * max_dist + 1) ^ 2

    nbhd0 <- make_nbhd(i = seq_len(nrow(env01[[2]])), sz = buffersize, 
                    r = env01[[1]], rdf = env01[[2]]) 

    message("Fitting model parameters")
    done <- list.files("data/output/simulations", pattern = "par_out_")
    done <- gsub("par_out_", "", done) %>%
            gsub(".rds", "", .) %>%
            as.numeric()
    todo <- setdiff(1:sim_n, done)
    message(paste0("Fitting ", length(todo), " individuals"))
    # fit <- do.call(rbind, lapply(todo, function(i) {
    foreach(i = todo, .combine = rbind) %dopar% {
        message(paste0("Fitting individual #: ", i, " / ", length(todo)))
        
        env01 <- list(terra::rast("data/output/simulations/env01.tif"),
                    readRDS("data/output/simulations/env01.rds"))

        current_jag <- i # for use in loglike_fun
        traject <- jag_traject_cells[[i]]
        prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
                rdf = env01[[2]])
        env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
            # Normalizing desired environmental variables for extended neighborhood
        env1 <- env1[nbhd_index]
        # Make indexing consistent with env
        names(env1) <- seq_len(length(nbhd_index))
        
        message("Fitting parameters for model 1: path-dependent kernel")
        objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
        ll <- log_likelihood(par0, objects1)
        message(paste0("Saving log-likelihood for model 1: ", i))
        saveRDS(ll, file = paste0("data/output/simulations/ll_fit1", i, ".rds"))
        par_out1 <- optim(par0, log_likelihood, objects = objects1)

        message("Fitting parameters for model 2: traditional SSF")
        obs_points <- as.data.frame(jag_traject[[i]])
        names(obs_points) <- c("x", "y")
        tr <- amt::steps(amt::make_track(obs_points, x, y))
        sl_emp <- as.vector(na.exclude(tr$sl_))
        ta_emp <- as.vector(na.exclude(tr$ta_))
        mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, max_dist = max_dist,
                                scale = 1)

        objects2 <- list(env1, max_dist, mk, obs)
        ll <- log_likelihood0(par0, objects2)
        message(paste0("Saving log-likelihood for model 2: ", i))
        saveRDS(ll, file = paste0("data/output/simulations/ll_fit2", i, ".rds"))
        par_out2 <- optim(par0, log_likelihood0, objects = objects2)
        
        message(paste0("COMPLETED path #: ", i, " / ", sim_n))
        saveRDS(c(par_out1$par, par_out2$par), 
                file = paste0("data/output/simulations/par_out_", i, ".rds"))
    }
    # ))
}
