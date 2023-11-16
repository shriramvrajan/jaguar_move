source("R/00_functions.R")

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("K", "RW", "RWM", "trad1", "trad2")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods

# excl <- which(res$nmove < 300)
res <- res[1:54, ]
# res <- res[-excl, ]


r2 <- res[, c("aic_RWM", "aic_trad1")]
r2 <- r2[order(r2$aic_RWM), ]
r2$ind <- seq_len(nrow(r2))
r2 <- gather(r2, "model", "aic", -ind)
bad <- r2$ind[which(is.na(r2$aic))]
r2 <- r2[-bad, ]


p <- ggplot(r2, aes(ind, aic, fill = model)) +
     geom_bar(stat = "identity", position = "dodge")
plot(p)

res <- res[order(res$aic_RWM), ]
plot(res$aic_RWM, pch = 19)
points(res$aic_trad1, pch = 19, col = "red")
abline(0, 1)

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