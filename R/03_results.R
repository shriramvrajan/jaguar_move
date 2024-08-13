source("R/00_functions.R")

# Results ----------------------------------------------------------------------
jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("o_trad", "o_rwm", "o_newfunc", "s2", "s6", "s7")
res <- results_table(s)[[1]][1:10, ]  # [[2]] = parameter values
head(res)
param <- results_table(s)[[2]]

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

m1 <- readRDS("data/output/debug_trad/debug_4.rds")
m2 <- readRDS("data/output/debug_pp/debug_4.rds")
