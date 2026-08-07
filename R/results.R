rm(list = ls())
source("R/functions.R")
source("R/classes.R")
Rcpp::sourceCpp("R/propagate.cpp")

## Simulations =================================================================

# Jump vs grain study
res <- readRDS("simulations/jumpgrain_20260806_233834.rds")$summary


## Empirical ===================================================================

### All filtering decisions must be explicit here
res_process <- function(r, m = 1) {
  res0 <- r$res_table
  if (anyNA(res0$pp_aic)) res0 <- res0[which(!is.na(res0$pp_aic)), ]
  if (anyNA(res0$ss_aic)) res0 <- res0[which(!is.na(res0$ss_aic)), ]
  if (any(res0$pp_conv != 0)) res0 <- res0[which(res0$pp_conv == 0), ]
  res <- merge(res0, jag_meta, by = c("ID", "biome"))
  nmoves <- data.frame(ID = as.numeric(jag_meta$ID),
                       nmove = nmove_by_multiple(m))
  res <- merge(res, nmoves, by = "ID")
  res <- res[res$nmove.y >= 30, ]  # nmove.x is total moves from jag_meta
  stopifnot(!anyNA(res$ID))       
  print(paste(nrow(res), "individuals"))
  return(res)
}

res <- results_set$new(r_ss = "data/output/emp_ss_1o_m1.rds", 
                       r_pp = "data/output/emp_pp_1o_m1.rds", env_type = "1o") %>%
                        res_process(., m = 1)

## Aggregate analysis ----------------------------------------------------------

### AIC
summ <- data.frame(id = res$ID, ss_ll = res$ss_ll, pp_ll = res$pp_ll,
                  diff = (res$ss_aic - res$pp_aic), 
                  conv = res$pp_conv,
                  nmove = res$nmove.y,
                  biome = res$biome,
                  meand = res$meandist,
                  jump = res$pp_njump,
                  sched = as.numeric(gsub("[^0-9.]", "", res$Planned.Schedule)))
# summ <- summ[-which(summ$conv == 52), ]
summ$diff[summ$diff < -4] <- -4
summ <- summ[order(summ$diff), ]
summ
print(length(which(summ$diff > 2)) / nrow(summ))
print(length(which(summ$ss_ll > summ$pp_ll)) / nrow(summ))
hist(summ$ss_ll - summ$pp_ll, 40, col = rgb(0, 0, 1), border = NA)
abline(v = 0, lwd = 3, col = "red")


### Holdout analysis -----------------------------------------------------------
h <- readRDS("data/output/holdout_1o_m1_f0.6_2026-08-05.rds")
hold <- ggplot(h, aes(x = nll_ss, y = nll_pp)) +
    geom_point(aes(color = "blue")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) 
plot(hold)
hh <- merge(res[, c("ID","nmove.y")], h, by = "ID")
