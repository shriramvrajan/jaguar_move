rm(list = ls())
library(pheatmap)
library(scales)
source("R/functions.R")
source("R/classes.R")
Rcpp::sourceCpp("R/propagate.cpp")

## Data and functions ==========================================================

r1 <- results_set$new(r_ss = "data/output/emp_ss_2026-07-10.rds", 
                      r_pp = "data/output/emp_pp_2026-07-13.rds", env_type = "1o")
res0 <- r1$res_table
res0 <- res0[which(!is.na(res0$pp_aic)), ]
res <- merge(res0, jag_meta, by = c("ID", "biome"))

## ouuuaoo
jdt <- lapply(jag_meta$ID, function(i) jaguar$new(i)$track$dt)
names(jdt) <- jag_meta$ID
par(mfrow = c(3, 3))
for (id in names(jdt)) {
  hist(jdt[[id]], 100, col = "black", main = id, border = NA)
}

## Aggregate analysis ==========================================================

### AIC ------------------------------------------------------------------------   
p_aic <- ggplot(res, aes(x = ss_aic, y = pp_aic)) +
    geom_point(aes(color = nmove)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) +
    labs(x = "Step selection AIC", y = "Path propagation AIC")
plot(p_aic)

diff <- res$ss_aic - res$pp_aic
print(sort(diff))
diff[diff < -4] <- -4
hist(diff, 100)
length(which(diff > 2))

### Deviance explained ---------------------------------------------------------
dev_exp <- function(i) {
  print(i)
  jag <- jaguar$new(id = i, results = res[res$ID == i, ])
  nulls <- jag$calculate_null_ll()
  de_dist <- 1 - nulls$ll_dist / nulls$ll_unif
  de_ss <- 1 - res$ss_ll[res$ID == i] / nulls$ll_unif
  de_pp <- 1 - res$pp_ll[res$ID == i] / nulls$ll_unif
  return(data.frame(id = i, de_dist = de_dist, 
                      de_ss = de_ss, de_pp = de_pp))
}
# k_fit <- data.frame(
#   id = sort(res$ID),
#   k = sapply(sort(res$ID), function(i) {
#     print(i)
#     jag <- jaguar$new(id = i, results = res[res$ID == i, ])
#     nulls <- jag$calculate_null_ll()
#     return(nulls$k)
#   }))
# saveRDS(k_fit, "data/output/k_fitted_start_values.rds")
dev <- do.call(rbind, lapply(res$ID, dev_exp))
dev_plot <- data.frame(
  ID = res$ID,
  val = c(dev$de_dist, dev$de_ss - dev$de_dist, dev$de_pp - dev$de_ss),
  col = rep(c("gray", "red", "blue"), each = nrow(dev))
)

dev_plot$col <- factor(dev_plot$col, levels = c("blue", "red", "gray"))

ggplot(dev_plot, aes(x = ID,  y = val, fill = col)) +
  geom_col() + 
  scale_fill_identity() +
  ylim(0, 1)

### Holdout analysis -----------------------------------------------------------
h <- readRDS("data/output/holdout_1o_2026-07-09.rds")
h$biome <- res$biome.x
h$nmove <- res$nmove
ind <- is.na(h$nll_ss) | is.na(h$nll_pp)
h <- h[-which(ind), ]
hold <- ggplot(h, aes(x = nll_ss / nmove, y = nll_pp / nmove)) +
    geom_point(aes(color = as.factor(biome))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) 
plot(hold)

plot(res$nmove, h$nll_ss - h$nll_pp)
cor(res$nmove, h$nll_ss - h$nll_pp, use = "pairwise.complete.obs")


## Individual analysis =========================================================

ll_compare <- function(id, env_type = "2o") {
    j_i <- jaguar$new(id = id, results = res[res$ID == id, ])
    track <- j_i$get_track()[, c("longitude", "latitude")]
    cells <- j_i$track_cells
    land <- j_i$get_landscape()
    sl_emp   <- na.exclude(j_i$get_track()$sl)
    max_dist <- ceiling(1.1 * max(sl_emp) / 1000)
    # extended neighborhood
    enbhd <- unique(as.numeric(na.exclude(make_nbhd(brdf, j_i$track_cells, 
                                                  sz = max_dist))))

    n_env <- if (env_type == "2o") 13 else if (env_type == "1o") 7 else stop()
    ss_par <- as.numeric(j_i$results[, paste0("ss_par", seq_len(n_env))])
    pp_par <- as.numeric(j_i$results[, paste0("pp_par", seq_len(n_env))])

    ss_phi <- pointwise_env(ss_par, enbhd, brdf, scale_from = enbhd) #env_type defaults to "2o"
    pp_phi <- pointwise_env(pp_par, enbhd, brdf, scale_from = enbhd)

    ss_ras <- trim(to_raster(ss_phi, enbhd, brazil_ras[[1]]))
    pp_ras <- trim(to_raster(pp_phi, enbhd, brazil_ras[[1]]))

    sim <- j_i$calculate_ll(env_type)
    non_outliers <- as.numeric(names(sim$ss$p_obs)) 
    delta <- abs(sim$ss$ll_obs) - abs(sim$pp$ll_obs)

    range_d <- range(delta, na.rm = TRUE)
    midpoint <- (0 - range_d[1]) / (range_d[2] - range_d[1])
    f_palette <- gradient_n_pal(colours = c("red", rgb(0, 0, 0, 0.4), "cyan"),
                                values = c(0, midpoint, 1))
    bins <- f_palette(scales::rescale(delta))

    par(mfrow = c(1, 2))
    terra::plot(ss_ras, main = "Step selection", col = gray.colors(100))
    points(track, col = "black", pch = 19, cex = 0.5)
    points(track[non_outliers, ], col = "blue", pch = 19, cex = 0.8)   
    terra::plot(pp_ras, main = "Path propagation", col = gray.colors(100))
    points(track, col = "black", pch = 19, cex = 0.5)
    points(track[non_outliers, ], col = "blue", pch = 19, cex = 0.8)   

    return(sim)
}

# ind <- c(12, 22, 50, 54, 99, 117)
# diag <- list()
# for (i in ind) {
#     print(i)
#     comp <- ll_compare(ind)
#     out <- setNames(list(t(comp$ss$p_surface), comp$pp$p_surface), c("ss", "pp"))
#     diag <- c(diag, out)
#     rm(comp)
#     rm(out)
#     gc()
# } 
# saveRDS(diag, "data/output/diagnostics.rds")

# id <- 62

bb <- ll_compare(62)

## Grain? ======================================================================

# Grain of the traversed landscape: the range of an exponential variogram fit
# to one covariate over the neighbourhood. Coordinates are raster row/col, so 
# range is in cells (1 cell = 1 km), same units as simulation autocorr. range.
grain_from_cells <- function(track_cells,  buffer = 10) {
  cells <- unique(as.numeric(na.exclude(make_nbhd(brdf, track_cells, sz = buffer))))
  cells_e <- scale(brdf[cells, 1:6])
  cells_e[is.na(cells_e)] <- 0
  pca_res <- prcomp(cells_e, center = TRUE, scale. = TRUE)
  z       <- as.numeric(pca_res$x[, 1])  # First principal component
  samp <- data.frame(z   = z,
                     row = brdf$row[cells],
                     col = brdf$col[cells])
  samp <- samp[stats::complete.cases(samp), ]
  tryCatch({
    sp::coordinates(samp) <- ~ col + row              # cell indices as planar coordinates
    vg  <- gstat::variogram(z ~ 1, samp)              # empirical semivariance vs lag
    fit <- gstat::fit.variogram(vg, gstat::vgm(model = "Exp"))
    fit$range[fit$model == "Exp"]                     # the e-folding scale in cells
  }, error = function(e) NA_real_)
}

grain_tbl <- do.call(rbind, lapply(res$ID, function(id) {
  print(id)
  jag <- jaguar$new(id)
  tr  <- jag$get_track()
  data.frame(
    ID          = id,
    grain       = grain_from_cells(jag$track_cells, buffer = 10),
    sl_cells    = median(na.exclude(tr$sl)) / 1000,   # displacement per fix, in cells
    dt_med      = median(na.exclude(tr$dt)),          # fix interval (min), mostly constant
    n_obs       = length(jag$track_cells)
  )
}))

outlier_ids <- c(12, 22, 50, 54, 99, 117)
pl <- merge(grain_tbl,
            res[, c("ID", "biome", "ss_aic", "pp_aic", "ss_ll", "pp_ll", "pp_njump")],
            by = "ID")
pl$dAIC    <- pl$ss_aic - pl$pp_aic                 # > 0: path propagation favoured
pl$dll_obs <- (pl$ss_ll - pl$pp_ll) / pl$n_obs      # per-observation ll gain (both are NLLs)
pl$ratio   <- pl$sl_cells / pl$grain        # grain-lengths crossed per fix
pl$outlier <- pl$ID %in% outlier_ids

# (A) The empirical Figure 3: individuals scattered on the prediction plane.
p1 <- ggplot(pl, aes(sl_cells, grain)) +
  geom_point(aes(fill = dAIC, shape = outlier), size = 3,
             colour = "grey20", stroke = 0.3) +
  scale_shape_manual(values = c("FALSE" = 21, "TRUE" = 24), guide = "none") +
  scale_fill_viridis_c(option = "magma", name = expression(Delta*AIC),
                       trans = pseudo_log_trans(base = 10)) +   # tames the huge outliers
  labs(x = "median displacement per fix (cells)",
       y = "landscape grain: variogram range (cells)") +
  theme_minimal(base_size = 12)
p1
# (B) The collapse to one axis: advantage against grain-lengths per fix.
p2 <- ggplot(pl, aes(ratio, dll_obs)) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 5),
              se = TRUE, colour = "grey40") +
  geom_point(aes(colour = biome, shape = outlier), size = 2.5) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
  scale_x_log10() +
  labs(x = "grain-lengths traversed per fix  (displacement / grain)",
       y = "per-observation ll gain  (SS - PP)") +
  theme_minimal(base_size = 12)
p2
