rm(list = ls())
source("R/functions.R")
source("R/classes.R")

batch_aic      <- FALSE
batch_holdout  <- FALSE
individual     <- TRUE

file_ss <- "data/output/empirical_ss_qcbs.rds"
file_pp <- "data/output/empirical_results_pp_2026-04-27.rds"
name <- paste0(gsub("data/output/empirical_|.rds", "", file_ss), "___",
              gsub("data/output/empirical_|.rds", "", file_pp))

r_ss <- readRDS(file_ss)
r_pp0 <- readRDS(file_pp)

pp_missing <- NA
pp_exists <- setdiff(seq_len(82), pp_missing)
## add NA entries to r_pp list for missing individuals
r_pp <- vector("list", length = length(r_ss))
r_pp[pp_exists] <- r_pp0
if (length(pp_missing) > 0) r_pp[pp_missing] <- NA

res <- results_table(r_ss = r_ss, r_pp = r_pp)
res <- cbind(res, jag_meta[, -c("ID", "biome")])
res <- res[-which(res$pp_conv != 0 | res$ss_conv != 0 | res$regular_moves <= 50), ]

# Empirical batch results ======================================================

if (batch_aic) {
  plot_pdf(nm = paste0("figs/aic_", name, ".pdf"), x = 6, y = 4)
  ggplot(res, aes(x = ss_aic, y = pp_aic, label = ID)) +
    geom_point(aes(col = biome)) +
    # geom_text(aes(label = ID, x = ss_aic + 500), size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Step selection AIC", y = "Path propagation AIC") 
  dev.off()
}

if (batch_holdout) {
  ll_holdout <- readRDS(paste0("data/output/holdout_", name, ".rds"))
  names(ll_holdout) <- c("ll_ss_hold", "ll_pp_hold", "ID")
  res <- merge(res, ll_holdout, by = "ID", all.x = TRUE)

  plot_pdf(nm = paste0("figs/holdout_ll_", name, ".pdf"), x = 6, y = 4)
  ggplot(res, aes(x = ll_ss_hold - ll_pp_hold, y = ss_aic - pp_aic, label = ID)) +
    geom_point(aes(col = biome)) +
    # geom_text(aes(label = ID, x = ll_ss + 0.1), size = 2) +
    # geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "deltaAIC", y = "deltaHoldout")
  dev.off()
}

# Empirical individual analysis ================================================

if (individual) {
  ids <- unlist(jag_id)[pp_exists]
  # ids <- c(20, 81)
  for (id in ids) {
    print(paste0("Analyzing individual ", id)) 
    if (file.exists(paste0("figs/empirical/", id, ".pdf"))) {
      message(paste0("Figure for individual ", id, " already exists. Skipping..."))
      next
    }
    jag_i <- individual_analysis$new(id = id, r_ss = r_ss, r_pp = r_pp)
    par_ss <- as.numeric(jag_i$results[3:11])
    par_pp <- as.numeric(jag_i$results[15:23])
    region <- get_local_region(jag_i$track_cells, brdf, buffer = 10)

    max_dist <- ceiling(1.1 * max(jag_i$track$sl, na.rm = TRUE) / 1000)
    scale_ref <- unique(as.vector(make_nbhd(rdf = brdf, i = jag_i$track_cells,
                                            sz = max_dist)))
    scale_ref <- scale_ref[!is.na(scale_ref)]

    rescale01 <- function(x) {
      r <- range(x, na.rm = TRUE)
      (x - r[1]) / (r[2] - r[1])
    }
    ## Pointwise surfaces
    phi_ss <- pointwise_env(par_ss, region, brdf, scale_ref) %>% rescale01()
    phi_pp <- pointwise_env(par_pp, region, brdf, scale_ref) %>% rescale01()
    phi_dist <- rescale01(phi_pp / sum(phi_pp, na.rm = TRUE))
    pp_model <- path_propagation_model$new()
    pi_pp <- pp_model$stationary_surface(par_pp, region, brdf, step_size = 1,
                                            scale_from = scale_ref) %>% rescale01()
    # log_ratio_pp <- log2(pi_pp / phi_pp)
    # log_ratio_pp[is.infinite(log_ratio_pp)] <- NA

    template <- brazil_ras[[1]]
    rast_phi_ss   <- to_raster(phi_ss, region, template)
    rast_phi_dist <- to_raster(phi_dist, region, template)
    rast_pi <- to_raster(pi_pp, region, template)
    # rast_log_ratio <- to_raster(log_ratio_pp, region, template)

    region_coords <- xyFromCell(brazil_ras, region)
    crop_ext <- ext(
      min(region_coords[, 1]), max(region_coords[, 1]),
      min(region_coords[, 2]), max(region_coords[, 2])
    )

    plot_pdf(nm = paste0("figs/empirical/", id, ".pdf"), x = 12, y = 4)
    par(mfrow = c(1, 2))
    terra::plot(crop(rast_phi_ss, crop_ext), main = "SSF: pointwise environment")
    points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    terra::plot(crop(rast_phi_dist, crop_ext), main = "PP: pointwise environment")
    points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    # terra::plot(crop(rast_pi,     crop_ext), main = "PP: stationary distribution")
    # points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    # terra::plot(crop(rast_log_ratio, crop_ext), main = id)
    # points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3))
    dev.off()
  }
}
