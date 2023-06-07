source("R/00_functions.R")

## Path generation parameters ==================================================

# x0 <- sample(100, 1); y0 <- sample(100, 1)
x0 <- 50
y0 <- 50
par0 <- c(-2, 2)           # env scaling parameter for states 1 & 2
tprob0 <- c(0.3, 0.6)      # transition probabilities 1→2, 2→1
step0 <- 3000              # number of steps to simulate
n_obs <- step0 / sim_interval  
# interval defined as global var in 00_functions.R

envsize <- 100 # size of landscape in cells
s1 <- 1        # strength of autocorrelation in cells
s2 <- 5
r1 <- 10       # range of autocorrelation in cells
r2 <- 30

## Landscape generation ========================================================

env01 <- gen_landscape(size = envsize, s = s1, r = r1)
env02 <- gen_landscape(size = envsize, s = s2, r = r2)
saveRDS(env01, "data/output/simulations/env01.RDS")
saveRDS(env02, "data/output/simulations/env02.RDS")

## Simulation ==================================================================

paths <- lapply(1:sim_n, function(i) {
          msg(paste0("Path #: ", i, " / ", sim_n))
          jag_path(x0, y0, step0, par = par0, neighb = buffersize, 
                   tprob = tprob0)
          })
saveRDS(paths, "data/output/simulations/paths.RDS")

## Fitting =====================================================================

# to_dest?
envdf <- env01[[2]]
env_index <- seq_len(nrow(envdf))
envdf <- cbind(envdf, env_index)
adj <- as.data.frame(raster::adjacent(env01[[1]], cells = env_index, 
                     include = TRUE, directions = 8, id = TRUE, sorted = TRUE))
to_dest <- do.call(rbind, lapply(env_index, function(i) {
    out <- adj$from[which(adj$to == i)]
    if (length(out) < 9) {
        out <- c(out, rep(NA, 9 - length(out)))
    }
    return(out)
}))

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
max_dist <- ceiling(max(unlist(dist)) * 2)
step_range <- (2 * max_dist + 1) ^ 2

nbhd0 <- make_nbhd(i = seq_len(nrow(env01[[2]])), sz = buffersize, r = env01[[1]], rdf = env01[[2]]) 

for (i in 1:sim_n) {
    current_jag <- i # for use in loglike_fun
    traject <- jag_traject_cells[[i]]
    input_prep(traject, max_dist, steps, nbhd0, r = env01[[1]], rdf = env01[[2]])
    
    env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
    env2 <- scales::rescale(env02[[2]]$sim1[nbhd_index], to = c(0, 1))

    # msg(paste0("Fitting path #: ", i, " / ", n))
    # par_optim <- rnorm(1)
    ll <- loglike_fun(par0[1])
    saveRDS(ll, file = paste0("data/output/simulations/ll", i, ".rds"))
    # msg(paste0("Done fitting path #: ", i, " / ", n))
}


## Testing =====================================================================

# par(mfrow = c(n, 2))
# for (i in seq_len(n)) {
#     p <- paths[[i]]
#     a1_t <- p$a1[which(p$state == 1)]
#     a1_f <- p$a2[which(p$state == 1)]
#     a2_t <- p$a2[which(p$state == 2)]
#     a2_f <- p$a1[which(p$state == 2)]

#     breaks1 <- (0:25 * 0.04) * max(c(a1_t, a1_f), na.rm = TRUE)
#     hist(a1_t, 30, col = rgb(0, 0, 1, 0.4), border = NA, breaks = breaks1)
#     abline(v = mean(a1_t, na.rm = TRUE), col = "blue")
#     hist(a1_f, 30, col = rgb(1, 0, 0, 0.4), add = TRUE, border = NA, 
#     breaks = breaks1)
#     abline(v = mean(a1_f, na.rm = TRUE), col = "red")

#     breaks2 <- (0:25 * 0.04) * max(c(a2_t, a2_f), na.rm = TRUE)
#     hist(a2_f, 30, col = rgb(1, 0, 0, 0.4), border = NA, breaks = breaks2)
#     abline(v = mean(a2_f, na.rm = TRUE), col = "red")
#     hist(a2_t, 30, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA,  
#          breaks = breaks2)
#     abline(v = mean(a2_t, na.rm = TRUE), col = "blue")
# }