rm(list = ls())
source("R/functions.R")
source("R/classes.R")
Rcpp::sourceCpp("R/propagate.cpp")

## Data and functions ==========================================================

jag_meta$regular_moves <- sapply(seq_len(nrow(jag_meta)), function(i) {
    print(i)
    jag_meta$nmove[i] - length(jaguar$new(jag_meta$ID[i], 
                                max_multiple = 2)$outliers)
  })

r1 <- results_set$new(r_ss = "data/output/emp_ss_1o_mm2.rds", 
                      r_pp = "data/output/emp_pp_2026-07-26.rds", env_type = "1o")


res0 <- r1$res_table
if (anyNA(res0$pp_aic)) res0 <- res0[which(!is.na(res0$pp_aic)), ]
if (any(res0$pp_conv != 0)) res0 <- res0[which(res0$pp_conv == 0), ]
res <- merge(res0, jag_meta, by = c("ID", "biome"))
res <- res[which(res$regular_moves >= 100), ]
print(paste(nrow(res), "individuals"))

## Aggregate analysis ==========================================================

### AIC ------------------------------------------------------------------------   
p_aic <- ggplot(res, aes(x = ss_aic, y = pp_aic)) +
    geom_point(aes(color = nmove)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) +
    labs(x = "Step selection AIC", y = "Path propagation AIC")
plot(p_aic)

summ <- data.frame(id = res$ID, ss_ll = res$ss_ll, pp_ll = res$pp_ll,
                  diff = (res$ss_aic - res$pp_aic) / res$regular_moves, 
                  conv = res$pp_conv)
# summ <- summ[-which(summ$conv == 52), ]
summ$diff[summ$diff < -4] <- -4
summ <- summ[order(summ$diff), ]
summ

### Deviance explained ---------------------------------------------------------


cat("SS worse than distance:", sum(dev$de_ss < dev$de_dist), "\n",
    "PP worse than SS:     ", sum(dev$de_pp < dev$de_ss),   "\n")

dev_plot <- data.frame(
  id   = rep(dev$id, 3),
  val  = c(dev$de_dist, dev$de_ss - dev$de_dist, dev$de_pp - dev$de_ss),
  part = factor(rep(c("distance", "+ environment", "+ path"), each = nrow(dev)),
                levels = c("+ path", "+ environment", "distance")))
dev_plot$id <- factor(dev_plot$id, levels = dev$id[order(dev$de_pp)])

ggplot(dev_plot, aes(x = id, y = val, fill = part)) +
  geom_col(width = 0.85) +
  scale_fill_manual(name = NULL,
    values = c(distance = "grey70", `+ environment` = "firebrick",
               `+ path` = "steelblue")) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "deviance explained") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, size = 6))
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

bb <- ll_compare(110)

id <- 117
jj <- jaguar$new(id, results = res[res$ID == id, ])
ppnew <- jj$results[, paste0("pp_par", seq_len(9))] %>% as.numeric
ppnew[9] <- -5
ll <- jj$calculate_ll("1o")
ll2 <- jj$calculate_ll(env_type = "1o", par_pp = ppnew)
plot(ll2$ss$ll_obs, ll2$pp$ll_obs)
abline(0, 1)

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
