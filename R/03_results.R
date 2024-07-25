source("R/00_functions.R")

plotstuff <- TRUE

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("holdRWM", "newfunc")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods

pp <- param[[2]]

plot1 <- function(i, p) {
    names <- c("footprint", "elevation", "slope", "forestcover", "dist2water", "dist2road")
    x <- seq(from = -1, to = 1, length.out = 100)
    par <- pp[i, c(2 * p, 2 * p + 1)]
    plot(x, par[1] * x + par[2] * x ^2, main = names[p], type = "l")
}

plotmain <- function(i) {
    par(mfrow = c(2, 3))
    for (p in 1:6) {
        plot1(i, p)
    }
}

# map_track and plotmain to compare