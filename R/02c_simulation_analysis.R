source("R/01b_simulations.R")

paths <- readRDS("data/output/simulations/paths.RDS")
probs <- load_if_exists(paste0("p", 1:sim_n, ".RDS"), 
                        dir = "data/output/simulations")
currents   <- load_if_exists(paste0("current", 1:sim_n, ".RDS"), 
                        dir = "data/output/simulations")

env01 <- rast("data/output/simulations/env01.tif")
env01 <- list(env01, raster_to_df(env01))
env02 <- rast("data/output/simulations/env02.tif")
env02 <- list(env02, raster_to_df(env02))

par(mfrow = c(1, 3))
res <- lapply(1:sim_n, function(i) {
    path <- paths[[i]]
    prob <- probs[[i]]
    curr1  <- currents[[i]][[1]]
    curr2  <- currents[[i]][[2]]
    side <- sqrt(dim(curr1)[1])

    state <- path$state[seq(1, nrow(path), sim_interval)]
    p     <- prob[nrow(prob), ]

    s <- 99
    ss <- 10
    r1 <- raster(nrows = side, ncols = side, vals = curr1[, s, ss])
    r2 <- raster(nrows = side, ncols = side, vals = curr2[, s, ss])
    terra::plot(r1, main = paste0("Kernel #", i))
    terra::plot(r2, main = paste0("Kernel #", i))
    hist(p[state == 1], 11, main = paste0("Path #", i), col = rgb(0, 0, 1, 0.5),
         xlim = c(0, max(p, na.rm = TRUE)), border = NA)
    hist(p[state == 2], 11, add = TRUE, col = rgb(1, 0, 0, 0.5), border = NA)
    abline(v = mean(p[state == 1], na.rm = TRUE), col = "blue")
    abline(v = mean(p[state == 2], na.rm = TRUE), col = "red")
    abline(v = 1 / (21 ^ 2))

    # plot_path(path)

    # x <- readline("Press [enter] to continue") 

    # a1_t <- path$a1[which(path$state == 1)]
    # a1_f <- path$a2[which(path$state == 1)]
    # a2_t <- path$a2[which(path$state == 2)]
    # a2_f <- path$a1[which(path$state == 2)]

    # breaks1 <- (0:25 * 0.04) * max(c(a1_t, a1_f), na.rm = TRUE)
    # hist(a1_t, 30, col = rgb(0, 0, 1, 0.4), border = NA, breaks = breaks1)
    # abline(v = mean(a1_t, na.rm = TRUE), col = "blue")
    # hist(a1_f, 30, col = rgb(1, 0, 0, 0.4), add = TRUE, border = NA, 
    # breaks = breaks1)
    # abline(v = mean(a1_f, na.rm = TRUE), col = "red")

    # breaks2 <- (0:25 * 0.04) * max(c(a2_t, a2_f), na.rm = TRUE)
    # hist(a2_f, 30, col = rgb(1, 0, 0, 0.4), border = NA, breaks = breaks2)
    # abline(v = mean(a2_f, na.rm = TRUE), col = "red")
    # hist(a2_t, 30, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA,  
    #      breaks = breaks2)
    # abline(v = mean(a2_t, na.rm = TRUE), col = "blue")
})
