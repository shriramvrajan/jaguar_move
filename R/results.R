rm(list = ls())
library(pheatmap)
source("R/functions.R")
source("R/classes.R")

r1 <- results_set$new(r_ss = "data/output/empirical_ss_qcbs.rds",
                      r_pp = "data/output/empirical_pp_qcbs.rds",
                      env_type = "1o")

# First-order results
r2 <- results_set$new(r_ss = "data/output/emp_ss_2026-05-22.rds",
                      r_pp = "data/output/emp_pp_2026-05-27.rds",
                      env_type = "2o")


diag1 <- r1$get_individual(id)$calculate_ll("1o")
diag2 <- r2$get_individual(id)$calculate_ll("2o")

res <- r2$res_table
# res <- res[!is.na(res$pp_conv) & res$pp_conv != 52, ]
p <- ggplot(res, aes(x = ss_aic, y = pp_aic)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) +
    labs(x = "Step selection AIC", y = "Path propagation AIC")
p


j <- jaguar$new(id = 13)

bb <- brdf[unique(j$track_cells), 1:6]
bb <- scale(bb[, 1:6],
                      center = colMeans(bb[, 1:6]),
                      scale  = apply(bb[, 1:6], 2, sd))
bb2 <- bb ^ 2
