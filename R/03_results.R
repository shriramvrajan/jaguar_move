source("R/00_functions.R")

plotstuff <- TRUE

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("holdRWM", "newfunc")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods


if (plotstuff) {
    res <- cbind(res, param$RWM)

    excl <- which(res$nmove < 100)
    res <- res[-excl, ]

    ggplot(res, aes(x = aic_holdRWM - aic_trad1, y = nmove, col = meandist)) + geom_point(size = 3)
}
