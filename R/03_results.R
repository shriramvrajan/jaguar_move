source("R/00_functions.R")

jagmeta_br <- jag_meta[ID %in% jag_id[[1]], ]
env <- brazil_ras

s <- c("null", "RWH", "RWM", "trad1")

res <- results_table(s)[[1]]  # [[2]] = parameter values
param <- results_table(s)[[2]]
# need to figure out what happened to the null likelihoods

res <- cbind(res, param$RWM)

excl <- which(res$nmove < 100)
res <- res[-excl, ]

scatter <- ggplot(res, aes(x = seq_len(nrow(res)), y = exp01(V7)))
scatter + geom_point(size = 3, aes(col = bio)) # + 
          geom_abline(intercept = 0, slope = 1, linetype = "dashed")

hist <- ggplot(res, aes(x = exp01(V7))) + geom_histogram(binwidth = 0.05)


ggplot(res, aes(x = weight, y = age, col = bio)) + geom_point(size = 3)

