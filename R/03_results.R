source("R/00_functions.R")

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

res <- results_table()[[1]]  # [[2]] = parameter values
param <- results_table()[[2]]
# need to figure out what happened to the null likelihoods

excl <- which(res$nmove < 300)
res <- res[-excl, ]


plot(res$aic_K, col = rgb(1, 1, 1, 0.4), pch = 19, ylim = c(0, 1.1e5))
points(res$aic_RW, col = rgb(0, 1, 0, 0.4), pch = 19)
points(res$aic_RWM, col = rgb(1, 0, 0, 0.4), pch = 19)
points(res$aic_trad2, col = rgb(0, 0, 1, 0.4), pch = 19)
points(res$aic_RWH, col = rgb(0.5, 0.5, 0, 0.4), pch = 19)

hist(res$aic_K, 10, col = rgb(0, 0, 0, 0.4), border = NA, xlim = c(0, 1.1e5))
hist(res$aic_RWM, 10, col = rgb(0, 1, 0, 0.4), add = TRUE, border = NA)



# only for trad
# res$p7 <- exp(res$p7) / (exp(res$p7) + 1)

# # hist(res$p7)
# plot(res$aic_rwh, res$aic_trad2, pch = 19)
# abline(0, 1)
# points(res$aic_rwh, res$aic_trad, pch = 19, col = "red")

# # hist(res$aic_rwh - res$aic_trad)


# # res <- res[-which(res$nmove < 100), ]

# message("Results compiled")

# if (FALSE) {
#    hist(aic_k, 100, col = rgb(1, 0, 0, 0.4), border = NA)
#     hist(aic_rw, 100, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA)
#     hist(aic_rwh, 100, col = rgb(0, 1, 0, 0.4), add = TRUE, border = NA)

#     plot(aic_k, pch = 19)
#     points(aic_rw, pch = 19, col = "blue")
#     points(aic_rwh, pch = 19, col = "green")
#     points(aic_h, pch = 19, col = "red")

#     barplot(aic_k, col = rgb(1, 0, 0, 0.4), border = NA)
#     barplot(aic_rw, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA)
#     barplot(aic_rwh, col = rgb(0, 1, 0, 0.4), add = TRUE, border = NA)
#     barplot(aic_h, col = rgb(0.5, 0.5, 0, 0.4), add = TRUE, border = NA)

# }