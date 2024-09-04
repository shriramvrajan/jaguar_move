source("R/00_functions.R")

# probably need to rewrite results_table function
s <- c("sim_ss", "sim_pp")

res <- results_table(s)

# jaguar
# hypothesis: model works better (aic_ss>aic_pp) if steps are longer
aicdiff <- res$aic_ss - res$aic_pp
distprop <- by(jag_move, jag_move$ID, function(x) {
    return(length(which(x$dist > 3000)) / length(x$dist))
}) %>% unlist() %>% as.numeric()
plot(distprop, )

jag_move$step_time <- as.POSIXct(jag_move$timestamp, format = "%m/%d/%y %H:%M")
jag_move <- jag_move[-which(ID == 25), ]
interval  <- by(jag_move, jag_move$ID, function(x) {
    return(difftime(x$step_time, lag(x$step_time), units = "hours"))
}) %>% unlist() %>% as.numeric()




# lets look at individuals first
i <- 52
lss <- paste0("data/output/sim_ss/out_", i, ".rds") %>% readRDS()
lpp <- paste0("data/output/out_", i, ".rds") %>% readRDS()
print(lss$out * 2 + 14)
print(lpp$out * 2 + 16)
l1 <- lss$attract
l2 <- t(lpp$array[, , which.max(rowSums(log(lpp$predictions), na.rm = T))])
plotl <- function(i) {
    plot(l1[i, ], type = "l", col = "red", ylim = c(0, 1))
    lines(l2[i, ], col = "blue")
}
plotl(330)

# now next step- compare AICs and look at how they relate to step lengths
results_table <- function(s) {

}




# Results ----------------------------------------------------------------------
jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

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
