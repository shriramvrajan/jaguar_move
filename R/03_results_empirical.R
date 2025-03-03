source("R/00_functions.R")

simnames <- c("ss", "ssh", "pp3", "pp3h")
res <- results_table(simnames)
# res0 <- res
# res <- res0[res0$nmove > 200, ]

ggplot(data = res, aes(x = aic_ss, y = aic_pp3)) +
    geom_point(mapping = aes(col = bio), size = 3) +
    geom_text(aes(label = ifelse(aic_ss < aic_pp3, jag_id, ""))) +
    geom_abline() 
hist(res$aic_ss - res$aic_pp3, breaks = 50)
abline(v = 0)

# Holdout set fitting ----------------------------------------------------------


# Old  -------------------------------------------------------------------------
## 6 jan 2025
# plot(res$aic_pp3 - res$aic_ss, res$nmove, pch = 19, col = "black", 
#      xlim = c(-500, 500))
# abline(v = 0)
# points(res$aic_pp3_2 - res$aic_ss, res$nmove, pch = 19, col = "red")
# points(res$aic_pp4 - res$aic_ss, res$nmove, pch = 19, col = "green")

## 16 dec 2024
# r0 <- res$aic_ss3 - res$aic_circ # -ve means L better, +ve R
# r1 <- res$aic_circpp - res$aic_circ
# r2 <- res$aic_circpp2 - res$aic_circ


# select <- sapply(seq_along(res$jag_id), function(x) {
#     options <- c(res$aic_circpp[x], res$aic_circpp2[x])
#     return(min(options))
# })

# r4 <- select - res$aic_circ
# r5 <- select - res$aic_ss3

# plot(select, pch = 19, col = rgb(1, 0, 0, 0.5))
# points(res$aic_circ, pch = 19, col = rgb(0, 1, 0, 0.5))

# # male vs female
# plotpdf(nm = "figs/malefemale.pdf", x = 6, y = 4)
# ggplot(data = res) +
#     geom_point(mapping = aes(x = aic_pp, y = aic_circpp, col = bio)) +
#     geom_abline()
# dev.off()

# plotpdf()
# plot(res$aic_ss3, res$aic_pp, pch = 19, cex = 0.7, xlab = "AIC for step selection",
#      ylab = "AIC for path propagation", col = "black")
# abline(0, 1, lty = 2)
# points(res$aic_ss3, res$aic_circpp, pch = 19, cex = 0.7, col = "red")
# # points(res$aic_ss3, res$aic_pp3, pch = 19, cex = 0.7, col = "#a03a15")
# dev.off()


# files <- paste0("data/output/o_trad/likelihood_", 1:82, ".rds")
# llt <- sapply(files, function(f) {
#     if(file.exists(f)) return(readRDS(f)) else return(NA)
# }) #trad empirical kernel
# names(llt) <- paste0("j", c(1:82))

# # jaguar
# # hypothesis: model works better (aic_ss>aic_pp) if steps are longer
# aicdiff <- res$aic_ss2 - res$aic_pp
# distprop <- by(jag_move, jag_move$ID, function(x) {
#     return(length(which(x$dist > 3000)) / length(x$dist))
# }) %>% unlist() %>% as.numeric()
# plot(distprop, aicdiff)
# plot(distprop, aicdiff, ylim = c(-500, 500))

# # jag_move$step_time <- as.POSIXct(jag_move$timestamp, format = "%m/%d/%y %H:%M")
# # jag_move <- jag_move[-which(ID == 25), ]
# # interval  <- by(jag_move, jag_move$ID, function(x) {
# #     return(difftime(x$step_time, lag(x$step_time), units = "hours"))
# # }) %>% unlist() %>% as.numeric()

# # lets look at individuals first
# i <- 52
# lss <- paste0("data/output/sim_ss/out_", i, ".rds") %>% readRDS()
# lpp <- paste0("data/output/out_", i, ".rds") %>% readRDS()
# print(lss$out * 2 + 14)
# print(lpp$out * 2 + 16)
# l1 <- lss$attract
# l2 <- t(lpp$array[, , which.max(rowSums(log(lpp$predictions), na.rm = T))])
# plotl <- function(i) {
#     plot(l1[i, ], type = "l", col = "red", ylim = c(0, 1))
#     lines(l2[i, ], col = "blue")
# }
# plotl(330)

# # now next step- compare AICs and look at how they relate to step lengths
# results_table <- function(s) {

# }




# # Results ----------------------------------------------------------------------
# jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
# env <- brazil_ras

# # Parameter plots --------------------------------------------------------------
# # need to figure out what happened to the null likelihoods
# # Plotting local landscape based on parameter values

# plot1 <- function(i, p, pp) {
#     names <- c("footprint", "elevation", "slope", "forestcover", "dist2water", "dist2road")
#     x <- seq(from = -1, to = 1, length.out = 100)
#     par <- pp[i, c(2 * p, 2 * p + 1)]
#     plot(x, par[1] * x + par[2] * x ^2, main = names[p], type = "l")
# }

# plotmain <- function(i, pp) {
#     par(mfrow = c(2, 3))
#     for (p in 1:6) {
#         plot1(i, p, pp)
#     }
# } 
# plotmain(1, param[[1]])

# # map_track and plotmain to compare

# # Debug ------------------------------------------------------------------------

# d1 <- readRDS("data/output/debug/debug_old_trad.rds")
# d2 <- readRDS("data/output/debug/debug_old_pp.rds")
# d3 <- readRDS("data/output/debug/debug_new_trad.rds")
# d4 <- readRDS("data/output/debug/debug_new_pp_fix.rds")

# l1 <- d1$like
# l2 <- d2$predictions[, c(-926, -927)][3, ]
# l3 <- d3$like
# l4 <- d4$predictions[, c(-926, -927)][3,]

# hist(l1, col = rgb(0, 0, 0, 0.3), border = NA, breaks = 50, xlim = c(0, 1))
# hist(l2, col = rgb(0.5, 0, 0.5, 0.3), border = NA, breaks = 50, add = TRUE)
# hist(l3, col = rgb(0, 0.5, 0.5, 0.3), border = NA, breaks = 30, xlim = c(0, 1))
# hist(l4, col = rgb(0.5, 0.5, 0, 0.3), border = NA, breaks = 30, add = TRUE)
