source("R/00_functions.R")

plotstuff <- FALSE

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("holdRWM", "holdtrad1")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods


if (plotstuff) {
    res <- cbind(res, param$RWM)

    excl <- which(res$nmove < 100)
    res <- res[-excl, ]

    ggplot(res, aes(x = aic_holdRWM - aic_holdtrad1, y = nmove, col = sex)) + geom_point(size = 3)
}
