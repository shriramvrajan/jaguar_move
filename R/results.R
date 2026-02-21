rm(list = ls())
source("R/functions.R")
source("R/classes.R")

# library(gridExtra)

file_ss <- "data/output/empirical_ss_qcbs.rds"
file_pp <- "data/output/empirical_pp_qcbs_rcpp.rds"
res <- results_table(file_ss = file_ss, file_pp = file_pp)
res <- cbind(res, jag_meta[, -c("ID", "biome")])

# get nmove without outliers - I don't need to be rerunning this every time!!
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
res <- res[-which(res$pp_conv != 0 | res$nmove1 <= 30), ]

batch_aic      <- TRUE
batch_holdout  <- FALSE
individual     <- FALSE

# Empirical batch results ======================================================

if (batch_aic) {
  plot_pdf(nm = "figs/aic2.pdf", x = 6, y = 4)
  ggplot(res, aes(x = ss_aic, y = pp_aic, label = ID)) +
    geom_point(aes(col = biome)) +
    geom_text(aes(label = ID, x = ss_aic + 500), size = 2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(x = "Step selection AIC", y = "Path propagation AIC") 
  dev.off()
}

if (batch_holdout) {
  ll_holdout <- sapply(seq_len(nrow(res)), function(i) {
    print(paste0("Holdout analysis for individual ", jag_id[i]))
    par_ss_i <- as.numeric(res[i, 3:11])
    par_pp_i <- as.numeric(res[i, 15:23])
    
  })
}

# Empirical individual analysis ================================================

if (individual) {
  # ids <- unlist(jag_id)
  ids <- c(20, 81)
  for (id in ids) {
    print(paste0("Analyzing individual ", id)) 
    if (file.exists(paste0("figs/empirical/", id, ".pdf"))) {
      message(paste0("Figure for individual ", id, " already exists. Skipping..."))
      next
    }
    jag_i <- individual_analysis$new(id = id, file_ss = file_ss, file_pp = file_pp)
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
    phi_dist <- phi_pp / sum(phi_pp, na.rm = TRUE) %>% rescale01()
    browser()
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
    par(mfrow = c(1, 3))
    terra::plot(crop(rast_phi_ss, crop_ext), main = "SSF: pointwise environment")
    points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    terra::plot(crop(rast_phi_dist, crop_ext), main = "PP: pointwise environment")
    points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    terra::plot(crop(rast_pi,     crop_ext), main = "PP: stationary distribution")
    points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.6, col = rgb(0, 0, 0, 0.3))
    # terra::plot(crop(rast_log_ratio, crop_ext), main = id)
    # points(cbind(jag_i$track$longitude, jag_i$track$latitude), pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3))
    dev.off()
  }
}


res2 <- readRDS("simulations/r1_obs_sweep_2026-02-13 13:20:41.184713.rds")$results
