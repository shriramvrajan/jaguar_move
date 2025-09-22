rm(list = ls())
source("R/functions.R")
source("R/classes.R")


# Empirical batch results ======================================================

res <- results_table(file_ss = "data/output/empirical_results_ss_2025-09-21.rds",
                     file_pp = "data/output/empirical_results_pp_2025-09-22.rds")

# Empirical individual analysis ================================================

jag_i <- individual_analysis$new(id = 115)
aa <- jag_i$compare_dispersal()
aa$step_selection$probs
aa$path_propagation$probs

# Simulation results ===========================================================

# frag <- readRDS("data/output/fragmentation_study_results_2025-08-08.rds")
# frag2 <- frag[-which(frag$ll_step == 0 | frag$ll_pp == 0), ]

# # plot(tapply(frag$ll_step - frag$ll_pp, frag$r1, mean))
# r1_values <- c(1, 2, 3, 4, 5, 6, 8, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100)
# plot(x = r1_values, y = tapply(frag2$ll_step - frag2$ll_pp, frag2$r1, mean), pch = 19)