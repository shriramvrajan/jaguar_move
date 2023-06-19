source("R/00_functions.R")

## Path generation parameters ==================================================

# x0 <- 50
# y0 <- 50
par0 <- c(1, 1)           # env scaling parameter for states 1 & 2
tprob0 <- c(0.5, 0.5)      # transition probabilities 1→2, 2→1
step0 <- 1000              # number of steps to simulate
n_obs <- step0 / sim_interval  
# interval defined as global var in 00_functions.R

envsize <- 200 # size of landscape in cells
s1 <- 6        # strength of autocorrelation 
s2 <- 1
r1 <- 5       # range of autocorrelation in cells
r2 <- 60

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
          msg(paste0("Path #: ", i, " / ", sim_n))
          # x0 <- sample(100, 1)
          # y0 <- sample(100, 1)
          x0 <- ceiling(envsize / 2)
          y0 <- ceiling(envsize / 2)
          
          jag_path(x0, y0, nstep = step0, par = par0, neighb = buffersize, 
                   tprob = tprob0)
          })
saveRDS(paths, "data/output/simulations/paths.RDS")

## Fitting =====================================================================

envdf <- env01[[2]]
env_index <- seq_len(nrow(envdf))
envdf <- cbind(envdf, env_index)

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
    msg(paste0("Testing path #: ", i, " / ", sim_n))
    current_jag <- i # for use in loglike_fun
    traject <- jag_traject_cells[[i]]
    input_prep(traject, max_dist, steps, nbhd0, r = env01[[1]], 
               rdf = env01[[2]])
    
    env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
    env2 <- scales::rescale(env02[[2]]$sim1[nbhd_index], to = c(0, 1))

    # par_optim <- rnorm(1)
    ll <- loglike_fun(par0[1])
    saveRDS(ll, file = paste0("data/output/simulations/ll", i, ".rds"))
}
