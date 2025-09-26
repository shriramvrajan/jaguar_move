rm(list = ls())
source("R/functions.R")
source("R/classes.R")

file_ss <- "data/output/empirical_results_ss_2025-09-24.rds"
file_pp <- "data/output/empirical_results_pp_2025-09-24.rds"

# Empirical batch results ======================================================

res <- results_table(file_ss = file_ss, file_pp = file_pp)

# Empirical individual analysis ================================================

jag_i <- individual_analysis$new(id = 114, file_ss = file_ss, file_pp = file_pp)
jag_i$compare_dispersal(step = 1, max_dist = 50)

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
