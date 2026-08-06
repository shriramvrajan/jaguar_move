rm(list = ls())
source("R/functions.R")
source("R/classes.R")
Rcpp::sourceCpp("R/propagate.cpp")

## Data and functions ==========================================================

nmove_by_multiple <- function(m) {
  sapply(seq_len(nrow(jag_meta)), function(i) {
      print(i)
      jag_meta$nmove[i] - 
        length(jaguar$new(jag_meta$ID[i], max_multiple = m)$outliers)
    })
}

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

## Aggregate analysis ==========================================================

### AIC ------------------------------------------------------------------------   

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


### Deviance explained ---------------------------------------------------------
dev <- do.call(rbind, lapply(res$ID, function(i) dev_exp(i, res = res)))
cat("SS worse than distance:", sum(dev$de_ss < dev$de_dist), "\n",
    "PP worse than SS:     ", sum(dev$de_pp < dev$de_ss),   "\n")


### Holdout analysis -----------------------------------------------------------
h <- readRDS("data/output/holdout_1o_m1_f0.6_2026-08-05.rds")
hold <- ggplot(h, aes(x = nll_ss, y = nll_pp)) +
    geom_point(aes(color = "blue")) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
    geom_text(aes(label = ID), vjust = -0.7, size = 3) 
plot(hold)
hh <- merge(res[, c("ID","nmove.y")], h, by = "ID")


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

# outlier_ids <- c(12, 22, 50, 54, 99, 117)
pl <- merge(grain_tbl,
            res[, c("ID", "biome", "ss_aic", "pp_aic", "ss_ll", "pp_ll", "pp_njump")],
            by = "ID")
pl$dAIC    <- pl$ss_aic - pl$pp_aic                 # > 0: path propagation favoured
pl$dll_obs <- (pl$ss_ll - pl$pp_ll) / pl$n_obs      # per-observation ll gain (both are NLLs)
pl$ratio   <- pl$sl_cells / pl$grain        # grain-lengths crossed per fix
pl$outlier <- pl$ID %in% outlier_ids

# (A) The empirical Figure 3: individuals scattered on the prediction plane.
p1 <- ggplot(pl, aes(sl_cells, grain)) +
  geom_point(aes(fill = dAIC), size = 3,
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
