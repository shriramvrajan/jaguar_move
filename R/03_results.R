source("R/00_functions.R")

# Results ----------------------------------------------------------------------
jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("sim_ss", "sim_pp")
# s <- c("o_trad", "o_rwm", "o_newfunc", "s2", "s6", "s7")
res <- results_table(s)[[1]]  # [[2]] = parameter values
head(res)
param <- results_table(s)[[2]]

plot(res$aic_sim_ss, res$aic_sim_pp, xlab = "AIC traditional", ylab = "AIC path propagation")
abline(0, 1)

# 

# plotting histograms for individuals 5 and 54
# our model is slightly worse for 5, better for 54p

plot2 <- function(i) {


    hist(pss, col = rgb(1, 0, 0, 0.3), border = NA, breaks = 50, xlim = c(0, 1))
    hist(ppp, col = rgb(0, 0, 1, 0.3), border = NA, breaks = 80, add = TRUE)
}










# Parameter plots --------------------------------------------------------------
# need to figure out what happened to the null likelihoods
# Plotting local landscape based on parameter values

plot1 <- function(i, p, pp) {
    names <- c("footprint", "elevation", "slope", "forestcover", "dist2water", "dist2road")
    x <- seq(from = -1, to = 1, length.out = 100)
    par <- pp[i, c(2 * p, 2 * p + 1)]
    plot(x, par[1] * x + par[2] * x ^2, main = names[p], type = "l")
}

plotmain <- function(i, pp) {
    par(mfrow = c(2, 3))
    for (p in 1:6) {
        plot1(i, p, pp)
    }
} 
plotmain(1, param[[1]])

# map_track and plotmain to compare

# Debug ------------------------------------------------------------------------

d1 <- readRDS("data/output/debug/debug_old_trad.rds")
d2 <- readRDS("data/output/debug/debug_old_pp.rds")
d3 <- readRDS("data/output/debug/debug_new_trad.rds")
d4 <- readRDS("data/output/debug/debug_new_pp_fix.rds")

l1 <- d1$like
l2 <- d2$predictions[, c(-926, -927)][3, ]
l3 <- d3$like
l4 <- d4$predictions[, c(-926, -927)][3,]

hist(l1, col = rgb(0, 0, 0, 0.3), border = NA, breaks = 50, xlim = c(0, 1))
hist(l2, col = rgb(0.5, 0, 0.5, 0.3), border = NA, breaks = 50, add = TRUE)
hist(l3, col = rgb(0, 0.5, 0.5, 0.3), border = NA, breaks = 30, xlim = c(0, 1))
hist(l4, col = rgb(0.5, 0.5, 0, 0.3), border = NA, breaks = 30, add = TRUE)
