source("R/00_functions.R")

## Parameters ==================================================================

### Landscape generation
envsize <- 100  # size of landscape in cells
s1 <- 6         # strength of autocorrelation 
s2 <- 1
r1 <- 5        # range of autocorrelation in cells
r2 <- 20

### Path generation
# x0 <- 50
# y0 <- 50

par0 <- c(2.5, -2)           # env scaling parameter for states 1 & 2
tprob0 <- c(0.1, 0.7)      # transition probabilities 1→2, 2→1

sim_interval <- 10   # GPS observations taken every n steps, for simulation
step0        <- 2000 # Number of steps to simulate
n_obs        <- step0 / sim_interval  
sim_steps    <- 10  # Number of steps to simulate
sim_n        <- 5   # Number of simulations to run

message(paste0("Simulation interval: ", sim_interval))
message(paste0("Observations per simulation: ", n_obs))
message(paste0("Simulation steps: ", sim_steps))
message(paste0("Number of simulations: ", sim_n))

## Landscape generation ========================================================

env01 <- gen_landscape(size = envsize, s = s1, r = r1)
env02 <- gen_landscape(size = envsize, s = s2, r = r2)

writeRaster(env01[[1]], "data/output/simulations/env01.tif", overwrite = TRUE)
writeRaster(env02[[1]], "data/output/simulations/env02.tif", overwrite = TRUE)

par(mfrow = c(1, 2))
terra::plot(env01[[1]])
terra::plot(env02[[1]])

## Simulation ==================================================================

paths <- lapply(1:sim_n, function(i) {
          message(paste0("Path #: ", i, " / ", sim_n))
          # x0 <- sample(100, 1)
          # y0 <- sample(100, 1)
          x0 <- ceiling(envsize / 2)
          y0 <- ceiling(envsize / 2)
          jag_path(x0, y0, step0, par = par0, neighb = buffersize, 
                   tprob = tprob0)
          })
saveRDS(paths, "data/output/simulations/paths.RDS")

plot_path(paths[[3]]) # red is state 1, blue is state 2

## Fitting =====================================================================

if (FALSE) {
    # to_dest? 
    envdf <- env01[[2]]
    env_index <- seq_len(nrow(envdf))
    envdf <- cbind(envdf, env_index)
    # adj <- as.data.frame(raster::adjacent(env01[[1]], cells = env_index, 
    #                      include = TRUE, directions = 8, id = TRUE, sorted = TRUE))
    # to_dest <- do.call(rbind, lapply(env_index, function(i) {
    #     out <- adj$from[which(adj$to == i)]
    #     if (length(out) < 9) {
    #         out <- c(out, rep(NA, 9 - length(out)))
    #     }
    #     return(out)
    # }))

    jag_traject <- lapply(paths, function(p) {
        out <- cbind(p$x, p$y, p$state)
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
    # max_dist <- ceiling(max(unlist(dist)) * 2)
    max_dist <- 10
    step_range <- (2 * max_dist + 1) ^ 2

    nbhd0 <- make_nbhd(i = seq_len(nrow(env01[[2]])), sz = buffersize, 
                    r = env01[[1]], rdf = env01[[2]]) 

    for (i in 1:sim_n) {
        message(paste0("Testing path #: ", i, " / ", sim_n))
        current_jag <- i # for use in loglike_fun
        traject <- jag_traject_cells[[i]]
        prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
                rdf = env01[[2]])
        
        env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
        env2 <- scales::rescale(env02[[2]]$sim1[nbhd_index], to = c(0, 1))

        # par_optim <- rnorm(1)
        ll <- log_likelihood(par0)
        saveRDS(ll, file = paste0("data/output/simulations/ll", i, ".rds"))
    }
}