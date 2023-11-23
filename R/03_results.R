source("R/00_functions.R")

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("null", "RW", "RWM", "trad1", "trad2")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods

excl <- which(res$nmove < 100)
res <- res[-excl, ]

res <- res[order(abs(res$aic_trad1 - res$aic_RWM), decreasing = TRUE), ]

scatter <- ggplot(res, aes(x = aic_RW - aic_RWM, y = log(totdist)))
scatter + geom_point(size = 3, aes(col = bio)) # + 
          geom_abline(intercept = 0, slope = 1, linetype = "dashed")


ggplot(res, aes(x = weight, y = age, col = bio)) + geom_point(size = 3)

plot(res$aic_trad2 - res$aic_RWM, res$nmove)
abline(h=1000)
