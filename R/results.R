rm(list = ls())
source("R_backup/functions.R")
source("R_backup/classes.R")

# library(gridExtra)

file_ss <- "data/output/empirical_results_ss_2025-11-26.rds"
file_pp <- "data/output/empirical_results_pp_2025-11-29.rds"

# Empirical batch results ======================================================

res <- results_table(file_ss = file_ss, file_pp = file_pp)
res <- cbind(res, jag_meta[, -c("ID", "biome")])

# get nmove without outliers
res$nmove1 <- sapply(as.vector(jag_id)$jag_id, function(i) {
  print(i)
  jag <- jaguar$new(as.numeric(i))
  track <- jag$get_track()
  n_obs <- length(jag$get_track_cells())
  dt_scaled <- track$dt[2:length(track$dt)] / median(na.exclude(track$dt))
  dt_discrete <- pmax(1, round(dt_scaled))
  outliers <- which(dt_discrete > 1)
  return(n_obs - length(outliers))
})
res2 <- res[-which(res$pp_conv != 0 | res$nmove1 <= 200), ]
res2$perf <- as.numeric(res2$ss_aic > res2$pp_aic)

plotpdf(nm = "figs/aic.pdf", x = 6, y = 4)
p1 <- ggplot(res, aes(x = ss_aic, y = pp_aic, label = ID)) +
  geom_point(aes(size = nmove1, col = as.numeric(ss_aic > pp_aic))) +
  # geom_text() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Step selection AIC", y = "Path propagation AIC") 
p1
dev.off()


res3 <- res2[-which(res2$nmove <= 200), ]

# plotpdf(nm = "figs/pars.pdf")
# plot(ggplot(res3, aes(x = ss_par9, y = pp_par9)) +
#   geom_point(aes(col = perf, size = nmove1)) + geom_abline())
# dev.off()

# Empirical individual analysis ================================================

# jag_i <- individual_analysis$new(id = 114, file_ss = file_ss, file_pp = file_pp)

# jag_i$compare_dispersal(step_size = 2, n_steps = 20)

# Simulation results ===========================================================

# frag <- readRDS("data/output/fragmentation_study_results_2025-08-08.rds")
# frag2 <- frag[-which(frag$ll_step == 0 | frag$ll_pp == 0), ]

# # plot(tapply(frag$ll_step - frag$ll_pp, frag$r1, mean))
# plot(x = r1_values, y = tapply(frag2$ll_step - frag2$ll_pp, frag2$r1, mean), pch = 19)

# max_dist <- 1000
# plotpdf(nm = paste0("figs/indiv_test_", Sys.Date(), ".pdf"), x = 8, y = 4)
# par(mfrow = c(1, 3))
# ssd <- jag_i$ss_disp$probs
# ssd[2002001] <- 0
# ssd <- ssd / max(ssd, na.rm = TRUE)
# rast_ss <- rast(matrix(ssd, 
#                 nrow = 2 * max_dist + 1, ncol = 2 * max_dist + 1))
# terra::plot(rast_ss, main = "Step selection")
# ppd <- jag_i$pp_disp$probs
# ppd[2002001] <- 0
# ppd <- ppd / max(ppd, na.rm = TRUE)
# rast_pp <- rast(matrix(ppd, 
#                 nrow = 2 * max_dist + 1, ncol = 2 * max_dist + 1))
# terra::plot(rast_pp, main = "Path propagation")
# diff <- ppd - ssd
# rast_diff <- rast(matrix(diff, 
#                 nrow = 2 * max_dist + 1, ncol = 2 * max_dist + 1))
# terra::plot(rast_diff, main = "Difference (PP - SS)")
# dev.off()
