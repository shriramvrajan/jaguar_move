rerun_sim <- FALSE

ifelse(rerun_sim, source("R/01b_simulations.R"), source("R/00_functions.R"))

paths <- readRDS("data/output/simulations/paths.RDS")

env01 <- rast("data/output/simulations/env01.tif")
env01 <- list(env01, raster_to_df(env01))

fit <- as.data.frame(readRDS("data/output/simulations/fit.RDS"))
names(fit) <- c("a1", "b1", "a2", "b2")

par(mfrow = c(1, 2))
plot(fit$a1, fit$b1)
abline(h = 0)
abline(v = 0)
plot(fit$a2, fit$b2)
abline(h = 0)
abline(v = 0)

sim_n <- nrow(fit)
ll1 <- -unlist(load_if_exists(paste0("ll_fit1", 1:sim_n, ".rds"), 
                             dir = "data/output/simulations"))
ll2 <- -unlist(load_if_exists(paste0("ll_fit2", 1:sim_n, ".rds"), 
                             dir = "data/output/simulations"))
aic1 <- 2 * 2 - 2 * ll1
aic2 <- 2 * 2 - 2 * ll2


par(mfrow = c(1, 2))
plot(aic1, aic2)
abline(0, 1)
hist(aic1 - aic2, breaks = 20)
abline(v = 0, lty = 2)
