source("R/00_functions.R")

## Landscape generation ========================================================

env1 <- gen_landscape(size = 100, s = 1, r = 10)
env2 <- gen_landscape(size = 100, s = 5, r = 30)

## Path generation parameters ==================================================

# x0 <- sample(100, 1); y0 <- sample(100, 1)
x0 <- 50
y0 <- 50
par0 <- c(-2, 2)
tprob0 <- c(0.1, 0.6)
neighb0 <- 1
step0 <- 1000
n <- 10 # Number of paths to simulate
# breaks0 <- c(0:25 * 0.04) # histogram

## Simulation ==================================================================

paths <- lapply(1:n, function(i) {
          msg(paste0("Path #: ", i, " / ", n))
          jag_path(x0, y0, step0, par = par0, neighb = neighb0, tprob = tprob0)
          })

## Fitting =====================================================================

step_range <- neighb0
n_obs <- step0
steps <- 25
interval <- 5 # GPS observations taken every n steps
buffersize <- 1 # How far does jaguar move in 1 time step

# to_dest?
envdf <- env1[[2]]
env_index <- seq_len(nrow(envdf))
envdf <- cbind(envdf, env_index)
adj <- as.data.frame(raster::adjacent(env1[[1]], cells = env_index, 
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
    ind <- seq(1, nrow(out), interval)
    out <- out[ind, ]
    return(out)
})

jag_traject_cells <- lapply(jag_traject, function(tr) {
    raster::cellFromXY(env1[[1]], tr[, 1:2])
})

dist <- lapply(jag_traject, function(tr) {
    out <- c(0, sqrt(diff(tr[, 1])^2 + diff(tr[, 2])^2))
    return(out)
})
max_dist <- ceiling(max(unlist(dist)) * 2)
step_range <- (2 * max_dist + 1) ^ 2

nbhd0 <- make_nbhd(i = seq_len(nrow(env1[[2]])), sz = buffersize, r = env1[[1]], rdf = env1[[2]]) 

for (i in 1:n) {
    # i <- 1
    traject <- jag_traject_cells[[i]]
    input_prep(traject, max_dist, steps, nbhd0, r = env1[[1]], rdf = env1[[2]])
    
    env1 <- scales::rescale(env1[[2]]$sim1[nbhd_index], to = c(0, 1))
    env2 <- scales::rescale(env2[[2]]$sim1[nbhd_index], to = c(0, 1))
    par_optim <- rnorm(1)
    o <- optim(par_optim, loglike_fun)
    saveRDS(o, file = paste0("data/output/simulation_optim_", i, ".rds"))
}


## Testing =====================================================================

# env11 <- exp(env1[[1]] * par0[1]) 
# env11 <- env11 / cellStats(env11, sum)
# env22 <- exp(env2[[1]] * par0[2])
# env22 <- env22 / cellStats(env22, sum)

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