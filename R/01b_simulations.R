source("R/00_functions.R")

## Landscape generation ========================================================

env1 <- gen_landscape(size = 100, s = 0.5, r = 30)
env2 <- gen_landscape(size = 100, s = 0.5, r = 30)
# env2 <- gen_landscape(size = 100, s = 0.033, r = 20)


## Path generation parameters ==================================================

# x0 <- sample(100, 1); y0 <- sample(100, 1)
x0 <- 50
y0 <- 50
par0 <- c(2, 2)
tprob0 <- c(0.5, 0.5)
neighb0 <- 6
step0 <- 1000
n <- 3 # Number of paths to simulate
# breaks0 <- c(0:25 * 0.04) # histogram

## Simulation ==================================================================

paths <- lapply(1:n, function(i) {
          msg(paste0("Path #: ", i, " / ", n))
          jag_path(x0, y0, step0, par = par0, neighb = neighb0, tprob = tprob0)
          })

## Testing =====================================================================

env11 <- exp(env1[[1]] * par0[1]) 
env11 <- env11 / cellStats(env11, sum)
env22 <- exp(env2[[1]] * par0[2])
env22 <- env22 / cellStats(env22, sum)

par(mfrow = c(n, 2))
for (i in seq_len(n)) {
    p <- paths[[i]]
    val1_t <- env11[p$cell][which(p$state == 1)]
    val1_f <- env11[p$cell][which(p$state == 2)]
    val2_t <- env22[p$cell][which(p$state == 2)]
    val2_f <- env22[p$cell][which(p$state == 1)]

    breaks1 <- (0:25 * 0.04) * max(c(val1_t, val1_f), na.rm = TRUE)
    hist(val1_t, 30, col = rgb(0, 0, 1, 0.4), border = NA, breaks = breaks1)
    abline(v = mean(val1_t, na.rm = TRUE), col = "blue")
    hist(val1_f, 30, col = rgb(1, 0, 0, 0.4), add = TRUE, border = NA, 
    breaks = breaks1)
    abline(v = mean(val1_f, na.rm = TRUE), col = "red")

    breaks2 <- (0:25 * 0.04) * max(c(val2_t, val2_f), na.rm = TRUE)
    hist(val2_f, 30, col = rgb(1, 0, 0, 0.4), border = NA, breaks = breaks2)
    abline(v = mean(val2_f, na.rm = TRUE), col = "red")
    hist(val2_t, 30, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA,  
         breaks = breaks2)
    abline(v = mean(val2_t, na.rm = TRUE), col = "blue")
}



## Fitting =====================================================================
if (FALSE) { 
step_range <- neighb0
n_obs <- step0
steps <- 25

# to_dest?
envdf <- env1[[2]]
env_index <- seq_len(nrow(envdf))
envdf <- cbind(envdf, env_index)
adj <- as.data.frame(raster::adjacent(env1[[1]], cells = env_index, include = T,
                                      directions = 8, id = T, sorted = T))
to_dest <- do.call(rbind, lapply(env_index, function(i) {
    out <- adj$from[which(adj$to == i)]
    if (length(out) < 9) {
        out <- c(out, rep(NA, 9 - length(out)))
    }
    return(out)
}))

# obs?
# Building observed data to test against
nbhd0 <- make_nbhd(r = env1[[1]], rdf = env1[[2]], i = 4950,
                   sz = neighb0)
index_mat <- matrix(
    data = seq_len(length(env_index)),
    nrow = (2 * neighb0 + 1)^2,
    ncol = n_obs
)
obs <- vector(length = ncol(index_mat) - 1)
for (y in 1:(ncol(index_mat) - 1)) {                                         # 13s
    test <- which(nbhd_index == jag_traject_cells[y + 1])
    num <- which(index_mat[1, y] < test & test < index_mat[nrow(index_mat), y])
    obs[y] <- which(index_mat[, y] == test[num])
} 

## Plotting ====================================================================

# raster::plot(env1[[1]])
# raster::plot(exp(par1[1] * env1[[1]]))
# raster::plot(env2[[1]])
# raster::plot(exp(par1[2] * env2[[1]]))

    path <- get(paste0("path", i))
    plot_path(path, new = F, vgram = F,  col = c("red", "blue")[path$state])
    points(x0, y0, pch = 19, cex = 0.8)

}