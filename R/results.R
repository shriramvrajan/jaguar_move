rm(list = ls())
source("R/functions.R")
source("R/classes.R")
Rcpp::sourceCpp("R/propagate.cpp")

## Simulations =================================================================

# Jump vs grain study
# res <- readRDS("simulations/jumpgrain_20260806_233834.rds")$summary
# p  <- ggplot(r, aes(s, ll_diff_per_obs, colour = factor(p))) +
#       geom_line(aes(group = interaction(p, gen_step)), alpha = 1) +
#       geom_point(aes(shape = factor(gen_step)), size = 4) +
#       scale_x_log10() +
#       labs(x = "substep length / autocorrelation range  (s)",
#             y = "mean log(L) advantage per obs",
#             colour = "substeps per fix (p)", shape = "jump (cells/tick)") +
#       theme_minimal()
# plot(p)

## Empirical ===================================================================

### Keeping this function deliberately outside of results_set class
### All filtering decisions must be made explicit here
res_process <- function(r, m = 1, obs_interval = 0) {
  res0 <- r$res_table
  if (anyNA(res0$pp_aic)) res0 <- res0[which(!is.na(res0$pp_aic)), ]
  if (anyNA(res0$ss_aic)) res0 <- res0[which(!is.na(res0$ss_aic)), ]
  if (any(res0$pp_conv != 0)) res0 <- res0[which(res0$pp_conv == 0), ]

  res <- merge(res0, jag_meta, by = c("ID", "biome"))
  nmoves <- data.frame(ID = as.numeric(jag_meta$ID),
                       nmove = nmove_by_multiple(m, obs_interval))
  res <- merge(res, nmoves, by = "ID")
  res <- res[res$nmove.y >= 30, ]  # nmove.y is jag_meta count minus outliers
  stopifnot(!anyNA(res$ID))       
  print(paste(nrow(res), "individuals"))
  return(res)
}

res <- results_set$new(r_ss = "data/output/emp_ss_1o_m1.rds", 
                       r_pp = "data/output/emp_pp_1o_m1.rds", env_type = "1o") %>%
                       res_process(., m = 1, obs_interval = 0)

res2 <- results_set$new(r_ss = "data/output/emp_ss_1o_gao.rds",
                        r_pp = "data/output/emp_pp_1o_m1_o0_2026-08-21.rds", env_type = "1o") %>%
                        res_process(., m = 1, obs_interval = 0)

ind1 <- which(res$ID %in% intersect(res$ID, res2$ID))
ind2 <- which(res2$ID %in% intersect(res$ID, res2$ID))
df <- data.frame(id = intersect(res$ID, res2$ID),
      ss = res$ss_ll[ind1], pp = res$pp_ll[ind1], nj = res$pp_njump[ind1],
      ssg = res2$ss_ll[ind2], ppg = res2$pp_ll[ind2], njg = res2$pp_njump[ind2])
de <- do.call(rbind, lapply(df$id, function(i) dev_exp(i, res = res2)))
df <- merge(df, de)
df <- df[order(df$de_pp - df$de_ss), ]

plot(df$nj, pch = 19, cex = 2, col = rgb(1, 0, 0, 0.5))
points(df$njg, pch = 19, cex = 2, col = rgb(0, 0, 1, 0.5))
cor(df$nj, df$njg)

ll_compare(117, res = res)
ll_compare(61, res = res2)

dev_exp(117, res = res)
dev_exp(117, res = res2)

j <- jaguar$new(41); o <- j$pp_objects(2)
p <- c(rnorm(7, sd = .1), -0.5, log(1), -20)
system.time(replicate(20, o$model$log_likelihood(p, o$objects, sim = FALSE, env_type = "1o")))

## Aggregate analysis ----------------------------------------------------------

### AIC
res <- res0
summ <- data.frame(id = res$ID, ss_ll = res$ss_ll, pp_ll = res$pp_ll,
                  diff = (res$ss_aic - res$pp_aic), 
                  conv = res$pp_conv,
                  nmove = res$nmove.y,
                  biome = res$biome,
                  meand = res$meandist,
                  jump = res$pp_njump,
                  sched = as.numeric(gsub("[^0-9.]", "", res$Planned.Schedule)))
summ$diff[summ$diff < -4] <- -4
summ <- summ[order(summ$diff), ]
summ
print(length(which(summ$diff > 2)) / nrow(summ))
print(length(which(summ$ss_ll > summ$pp_ll)) / nrow(summ))
hist(summ$ss_ll - summ$pp_ll, 40, col = rgb(0, 0, 1), border = NA)
abline(v = 0, lwd = 2, lty = 2, col = "gray")
abline(v = 2, lwd = 2, lty = 1, col = "red")

### Holdout analysis -----------------------------------------------------------

h <- readRDS("data/output/holdout_1o_m1_f0.6_2026-08-05.rds")
hold <- ggplot(h, aes(x = nll_ss, y = nll_pp)) +
    geom_point(aes(color = "blue")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) 
plot(hold)
hh <- merge(res[, c("ID","nmove.y")], h, by = "ID")
