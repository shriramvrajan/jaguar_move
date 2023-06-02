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
step0 <- 3000
n <- 3 # Number of paths to simulate
# breaks0 <- c(0:25 * 0.04) # histogram

## Simulation ==================================================================

paths <- lapply(1:n, function(i) {
          msg(paste0("Path #: ", i, " / ", n))
          jag_path(x0, y0, step0, par = par0, neighb = neighb0, tprob = tprob0)
          })

## Testing =====================================================================

# env11 <- exp(env1[[1]] * par0[1]) 
# env11 <- env11 / cellStats(env11, sum)
# env22 <- exp(env2[[1]] * par0[2])
# env22 <- env22 / cellStats(env22, sum)

par(mfrow = c(n, 2))
for (i in seq_len(n)) {
    p <- paths[[i]]
    a1_t <- p$a1[which(p$state == 1)]
    a1_f <- p$a2[which(p$state == 1)]
    a2_t <- p$a2[which(p$state == 2)]
    a2_f <- p$a1[which(p$state == 2)]

    breaks1 <- (0:25 * 0.04) * max(c(a1_t, a1_f), na.rm = TRUE)
    hist(a1_t, 30, col = rgb(0, 0, 1, 0.4), border = NA, breaks = breaks1)
    abline(v = mean(a1_t, na.rm = TRUE), col = "blue")
    hist(a1_f, 30, col = rgb(1, 0, 0, 0.4), add = TRUE, border = NA, 
    breaks = breaks1)
    abline(v = mean(a1_f, na.rm = TRUE), col = "red")

    breaks2 <- (0:25 * 0.04) * max(c(a2_t, a2_f), na.rm = TRUE)
    hist(a2_f, 30, col = rgb(1, 0, 0, 0.4), border = NA, breaks = breaks2)
    abline(v = mean(a2_f, na.rm = TRUE), col = "red")
    hist(a2_t, 30, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA,  
         breaks = breaks2)
    abline(v = mean(a2_t, na.rm = TRUE), col = "blue")
}

## Fitting =====================================================================

step_range <- neighb0
n_obs <- step0
steps <- 25
interval <- 10 # GPS observations taken every n steps

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
step_range <- ceiling(max(unlist(dist)) * 2)

# for (i in 1:n) {
    i <- 1
    input_prep(jag_traject_cells[[i]], step_range, steps)
    env1i <- norm_env(env1[[2]], nbhd_index)
    env2i <- norm_env(env2[[2]], nbhd_index)
    o <- optim(par0, loglike_fun)
# }
