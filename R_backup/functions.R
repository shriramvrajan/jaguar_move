## Basic functions and libraries 

# Libraries ====================================================================
library(terra)
library(tidyverse)
library(data.table)
library(gstat)
library(ctmm)
library(amt) 
library(lubridate)
library(metaRange)
library(pryr)
library(foreach)
library(doParallel)
library(viridis)
library(ggplot2)

# Data =========================================================================

# CRS definitions
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
epsg5880 <- "+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 
+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# Jaguar movement data, ID numbers, and metadata
jag_move <- readRDS("data/jag_data_BR.rds")
jag_id   <- readRDS("data/jag_list.rds")
njag     <- nrow(jag_id)
jag_meta0 <- data.table(read.csv("data/input/jaguars/jaguar_metadata2.csv"))
jag_meta  <- readRDS("data/jag_meta2.rds")

# RasterStack of environmental variables -- see 01_generate_data.R for details
brazil_ras <- rast("data/env_layers.grd")
# load("data/env_layers.RData") # same information as a dataframe ('brdf')
brdf <- readRDS("data/env_layers.rds")
# brdf <- brdf[, -1]  # Remove human footprint

# Brazil biomes shapefile
biome <- vect("data/input/Brazil_biomes/Brazil_biomes.shp")

# Empirical results
# e_results <- readRDS("data/output/oldresults_empirical.rds")

# Functions ====================================================================

# 0. Basic ---------------------------------------------------------------------

# Output message m both in console and in logfile f
message <- function(m, f = "data/output/run_log.txt") {
    m <- paste(format(Sys.time(), "%d.%m.%y  %R"), m)
    print(m)
    cat(m, file = f, append = TRUE, sep = "\n")
}

# Object size list
memory <- function() {
    for(x in ls(envir = .GlobalEnv)) {
        paste(object.size(get(x)) %>% format(units = "Mb"), x) %>% print()
    }
    print(paste("Total:", as.numeric(mem_used()) / 1e9))
}

# Output last n lines of logfile
log_tail <- function(n = 10, f = "data/output/run_log.txt") {
    tail(readLines(f), n)
}

# Set up parallel processing
parallel_setup <- function(n_cores = 4) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    message(paste0("Parallel processing set up with ", n_cores, " cores."))
}

# Turn raster into a 3-column dataframe with row, col, and value
raster_to_df <- function(r) {
    outdf <- as.data.frame(r)
    outdf <- cbind(outdf, rowColFromCell(r, seq_len(nrow(outdf))))
    names(outdf)[(ncol(outdf) - 1):ncol(outdf)] <- c("row", "col")
    return(outdf)
}

# Plot to pdf (for SSH)
plotpdf <- function(nm = "plot1.pdf", x = 4, y = 4) {
  pdf(nm, width = x, height = y)
}

# Return 9-cell neighborhood of cell in raster r
rast9 <- function(r, cell) {
    # cell: cell number in raster
    # r: raster object
    row <- rowFromCell(r, cell)
    col <- colFromCell(r, cell)
    expand <- expand.grid((row - 1):(row + 1), (col - 1):(col + 1))
    cell9 <- ext(r, cells = terra::cellFromRowCol(r, expand$Var1, expand$Var2))
    zoom(r, cell9)
}

plot_satellite <- function(bbox) {
    tryCatch({
        satellite_map <- ggmap::get_map(location = bbox, maptype = "satellite")
        p <- ggmap(satellite_map) + labs(title = "Dispersal kernel area") +
        theme_void()
        print(p)
    }, error = function(e) {
        message("Satellite map not available.")
        return(NULL)
    })
}

# Bresenham line algorithm
bresenham <- function(x0, y0, x1, y1, grid_size) {
  # Clip to grid boundaries
  x0 <- round(max(1, min(grid_size, x0)))
  y0 <- round(max(1, min(grid_size, y0)))
  x1 <- round(max(1, min(grid_size, x1)))
  y1 <- round(max(1, min(grid_size, y1)))
  dx <- abs(x1 - x0)
  dy <- abs(y1 - y0)
  err <- dx - dy # ?
  sx <- if(x0 < x1) 1 else -1
  sy <- if(y0 < y1) 1 else -1

  cells <- matrix(ncol = 2, nrow = 0)
  x <- x0
  y <- y0
  
  while (TRUE) {
      cells <- rbind(cells, c(x, y))
      if (x == x1 && y == y1) break
      e2 <- 2 * err
      if (e2 > -dy) {
          err <- err - dy
          x <- x + sx
      }
      if (e2 < dx) {
          err <- err + dx
          y <- y + sy
      }
      # Safety check to prevent infinite loops
      if (x < 1 || x > grid_size || y < 1 || y > grid_size) break
  }
  return(cells)
}

# 1. Movement model ============================================================

# Output neighborhood as matrix
# r:    Raster object
# i:    Index of cell in raster
# sz:   Size of neighborhood
make_nbhd <- function(r = brazil_ras, rdf = brdf, i, sz) {
  n_cells <- length(i)
  n_offsets <- (2 * sz + 1)^2

  offsets_row <- t(rep(-sz:sz, each = 2 * sz + 1))
  offsets_col <- t(rep(-sz:sz, times = 2 * sz + 1))
  rows_all <- rep(rdf$row[i], each = n_offsets) + rep(offsets_row, times = n_cells)
  cols_all <- rep(rdf$col[i], each = n_offsets) + rep(offsets_col, times = n_cells)

  cells_all <- terra::cellFromRowCol(r, rows_all, cols_all)
  mat <- matrix(cells_all, nrow = n_cells, ncol = n_offsets, byrow = TRUE)
  return(mat)
}

make_movement_kernel <- function(n = 10000, sl_emp, ta_emp, max_dist, 
                                 minimum = 0, scale = 1000) {
  
  avail <- data.frame(sl = sample(sl_emp, n, replace = TRUE),
                      ta = sample(ta_emp, n, replace = TRUE))
  avail$x <- avail$sl * cos(avail$ta) / scale # scales from meters to km
  avail$y <- avail$sl * sin(avail$ta) / scale 
  avail$xi <- sapply(avail$x, function(x) ifelse(x > 0, floor(x), ceiling(x)))
  avail$yi <- sapply(avail$y, function(y) ifelse(y > 0, floor(y), ceiling(y)))
  density <- tapply(avail$x, list(avail$yi, avail$xi), length)
  density <- reshape2::melt(density)
  names(density) <- c("x", "y", "n")
  size <- max_dist * 2 + 1
  out <- data.frame(x = rep(-max_dist:max_dist, times = size),
                    y = rep(-max_dist:max_dist, each = size))

  out <- merge(out, density, by = c("x", "y"), all.x = TRUE)
  out$n[is.na(out$n)] <- 0
  out$n <- out$n + minimum
  out$n <- out$n / sum(out$n)
  
  return(out$n)
}

# Normalize probabilities across neighbors of each cell
# v: vector of cell values
normalize_nbhd <- function(v, nbhd) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = TRUE)
}

# Attractiveness function for movement model
# env:      Environmental variable values
# par:      Parameters for functional form
# format:   TRUE to return as matrix, FALSE to return as vector
env_function <- function(env, par, nbhd = NULL, sim = FALSE) {

  par <- as.numeric(par)
  if (!sim) {
    # env <- scale(env)
    if (any(is.na(env))) env[which(is.na(env))] <- 0
    # First order with intercept
    attract <- 1 / (1 + exp(par[1] * env[, 1] + par[2] * env[, 2] +
                            par[3] * env[, 3] + par[4] * env[, 4] + 
                            par[5] * env[, 5] + par[6] * env[, 6] + par[7]))
  } else {
    attract <- 1 / (1 + exp(par[1] + par[2] * env + par[3] * env^2))
  }

  if (is.null(nbhd)) {
    return(attract)
  } else {
    attract <- matrix(attract[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
    attract <- attract / rowSums(attract, na.rm = TRUE)
    return(attract)  
  }
}

apply_kernel <- function(attract0, kernel, bg_rate = 0) {
  kernel <- kernel / sum(kernel, na.rm = T)
  na_mask <- is.na(attract0)
  attract0[na_mask] <- 0
  p <- attract0 * rep(kernel, each = nrow(attract0))
  p <- p + bg_rate - p * bg_rate
  p <- p / rowSums(p, na.rm = TRUE)
  p[na_mask] <- NA
  return(p)
}


# 2. Data exploration ----------------------------------------------------------

# Plot kernel for movement model
plot_circle_kernel <- function(max_dist, k_par, bg_rate) {
    kern <- calculate_dispersal_kernel(max_dispersal_dist = max_dist,
                                       kfun = function(x) {
                                        dexp(x, rate = k_par) + bg_rate
                                       })
    print(sum(kern))
    terra::plot(rast(kern))
    return(kern)
}

# Plot environmental layers for path 
# layer: layer number of raster
# bounds: bounding box of path
# path: path of jaguar
# grad: number of colors to use for path
# ras: raster stack of environmental layers 
plot_env <- function(layer, bounds, path, grad = 20, ras = brazil_ras) {
    terra::plot(crop(ras[[layer]], bounds), main = names(ras)[layer], 
                col = (gray.colors(50)))
    # pal <- colorRampPalette(c("#ff0000", "#1bbe29"))(grad)
    pal <- viridis(grad)
    path <- as.data.frame(path)
    path$col <- rep(pal, each = ceiling(nrow(path) / grad))[seq_len(nrow(path))]
    points(x = path$x, y = path$y, pch = 19, cex = 0.5, col = "red")
    # lines(x = path$x, y = path$y, col = rgb(0, 0, 0, 0.3), lwd = 0.3)
    # for (i in seq_len(nrow(path) - 1)){
    #     segments(path$x[i], path$y[i], path$x[i+1], path$y[i+1], 
    #              col = path$col[i])
    #     # segments(df$x[i], df$y[i], df$x[i+1], df$y[i+1], col=df$col[i])
    # }

}

# Plot time series (in track form) as 3D line plot (x, y, t)
# tr: track object
plot_xyt <- function(tr) {
    plot_ly(x = tr$x_, y = tr$y_, z = tr$t_, type = "scatter3d", 
        mode = "lines", line = list(width = 1))
}

# Map path of jaguar i
# i: ID number of jaguar
# grad: number of colors to use for path
# type: 1 is satellite map using ggmap, 2 is env layers
map_track <- function(i, grad = 20, type = 2, file = FALSE) {
    # type: 1 is satellite map using ggmap, 2 is env layers
    moves <- jag_move[ID == as.numeric(i)]

    # Get bounding box of all GPS measurements
    gps <- vect(moves[, 3:4], geom = c("longitude", "latitude"))
    bbox <- ext(gps)

    # Get total tracking period
    dates <- as.Date(sapply(moves$timestamp, function(dt) {
        strsplit(as.character(dt), " ")[[1]][1]
    }), format = "%m/%d/%y")
    period <- difftime(dates[length(dates)], dates[1])

    names(moves)[3:4] <- c("x", "y")
    path <- sp::SpatialPoints(coords = moves[, 3:4], sp::CRS("+init=epsg:4326"))
    # path2 <- as.data.frame(sp::spTransform(path, OpenStreetMap::osm()))

    if (file) plotpdf(nm = "track.pdf", x = 8, y = 4)
    if (type == 1) {
        bboxgg <- bbox[c(1, 3, 2, 4)]
        names(bboxgg) <- c("left", "bottom", "right", "top")
        map0 <- ggmap::get_map(location = bboxgg, maptype = "satellite")
        ggmap(map0) +
            geom_point(aes(x = x, y = y),
                data = as.data.frame(path),
                color = "red"
            )
    } else if (type == 2) {
        par(mfrow = c(2, 3))
        for (i in 1:6) {
            plot_env(layer = i, bounds = bbox, path = path, grad = grad)
        }
    }
    if (file) dev.off()
}

# Map home range of jaguar across environmental layers
# id: ID number of jaguar
# vgram: TRUE to plot variogram of path
map_homerange <- function(id, vgram = FALSE) {
  message(paste0("Mapping home range of jaguar ", id, "..."))
  path <- jag_move[ID == as.numeric(id)]
  path <- as.telemetry(path, timeformat = "auto")
  if (vgram == TRUE) {
    vpath <- variogram(path)
    plot(vpath)
  }
  guess <- ctmm.guess(path, interactive = FALSE)
  path_guess <- ctmm.fit(path, guess)
  kde <- akde(path, path_guess)
  par(mfrow = c(1, 1))
  plot(kde)
  plot(path, add = TRUE)
}

# Decompose time series (in telemetry form) into 28d chunks
# tel: telemetry object
# n:   which 28d chunk to plot
lunarize <- function(tel, n = 1) {
    lunar <- 28 * 24 * 60 * 60
    t0 <- tel$t[1]
    start <- t0 + (n - 1) * lunar
    end   <- t0 + n * lunar
    
    tel1 <- tel[tel$t %in% start:end, ]
    pgram <- periodogram(tel1)    
    vgram <- variogram(tel1)

    par(mfrow = c(1, 2))
    plot(pgram, max = TRUE, diagnostic = TRUE, pch = 19, cex = 0.8)
    plot(vgram)
}

# Calculate frequency of staying in the same cell
calc_move_freq <- function(dat) {
  # Calculate how often row n = row n+1
  dat <- as.data.frame(dat)
  names(dat) <- c("x", "y")
  dat$stay <- ifelse(dat$x == lag(dat$x) & 
                     dat$y == lag(dat$y), 1, 0)
  return(mean(dat$stay, na.rm = T))
}

# Plot movement kernel
plot_curve <- function(par, mu = 0, sd = 1, bounds = c(0, 10), values = FALSE) {
    # Plot functional form for 2op, mu/sd reverse normalization of env variable.
    x <- seq(bounds[1], bounds[2], 0.1)
    x0 <- (x - mu) / sd
    y <- 1 / (1 + exp(par[1] + par[2] * x0 + par[3] * x0^2))
    if (values) return(y)
    plot(x, y, type = "l", xlab = "Environmental variable", 
         ylab = "Attraction")
    abline(v = 0, lty = 2)
}


# 2b. Simulation ---------------------------------------------------------------

thicken_line <- function(line_cells, width, grid_size) {
    if (width <= 1) return(line_cells)
    
    # Start with original line
    all_cells <- line_cells
    
    # For width=2, add immediate 8-connected neighbors
    # For width=3, add neighbors within distance 1, etc.
    radius <- floor(width / 2)
    
    for (i in seq_len(nrow(line_cells))) {
        x <- line_cells[i, 1]
        y <- line_cells[i, 2]
        # Add all cells within radius
        for (dx in 0:radius) {
            for (dy in 0:radius) {
                new_x <- x + dx
                new_y <- y + dy
                # Check bounds
                if (new_x >= 1 && new_x <= grid_size && 
                    new_y >= 1 && new_y <= grid_size) {
                    all_cells <- rbind(all_cells, c(new_x, new_y))
                }
            }
        }
    }
    return(unique(all_cells))
}

# Generates a random field with a given correlation structure
# beta: beta parameter for gstat, s: sill, r: range, n: nugget
gen_landscape <- function(size = 100, beta = 1, s = 0.03, r = 10, n = 0,
                      b_density = 0, b_length = 20, b_value = 2, b_width = 2) {   

    # Autocorrelation model - base field
    xy <- expand.grid(1:size, 1:size)
    names(xy) <- c("x", "y")
    model <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = beta, 
                   model = vgm(psill = s, range = r, model = "Exp", nugget = n), 
                   nmax = 20)
    out <- predict(model, newdata = xy, nsim = 1)
    if (any(out < 0)) out[out < 0] <- 0
    out_mat <- matrix(out$sim1, nrow = size, ncol = size, byrow = TRUE)

    # Add barriers
    if (b_density > 0) {
      # Number of barriers per 10k cells
      n_barriers <- ceiling(b_density * (size ^ 2) / 10000)

      for (b in seq_len(n_barriers)) {
        x0 <- sample(1:size, 1)
        y0 <- sample(1:size, 1)
        angle <- runif(1, 0, 2 * pi)

        length <- rgamma(1, shape = 4, scale = b_length / 4)
        x1 <- x0 + length * cos(angle)
        y1 <- y0 + length * sin(angle)

        line_xy <- bresenham(x0, y0, x1, y1, size)
        thick_cells <- thicken_line(line_xy, b_width, size)

        for(i in seq_len(nrow(thick_cells))) {
          out_mat[thick_cells[i, 1], thick_cells[i, 2]] <- b_value
        }
      }
    }

    # Output as raster 
    out <- data.frame(x = rep(1:size, each = size), y = rep(1:size, times = size),
                      sim1 = as.vector(t(out_mat))) %>% rast
    plotpdf(nm = "figs/current_landscape.pdf")
    terra::plot(out)
    dev.off()
    return(out)
}

# Generate a jaguar path of n steps starting from (x0, y0) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x0, y0, n_step, par, neighb, env_raster) {
  
  # Pre-compute neighborhood lookup for entire raster
  rdf <- raster_to_df(env_raster)
  env_values <- rdf$sim1
  all_cells <- seq_len(nrow(rdf))
  nbhd_lookup <- make_nbhd(i =  all_cells, sz = neighb, r = env_raster, rdf = rdf)

  path <- matrix(NA, nrow = n_step, ncol = 4)
  current_cell <- cellFromRowCol(env_raster, x0, y0)
  path[1, ] <- c(x0, y0, current_cell, NA)

  for (i in 2:n_step) {
      nbhd_cells <- nbhd_lookup[current_cell, ]

      env_vals <- env_values[nbhd_cells]
      attract <- 1 / (1 + exp(par[1] + par[2] * env_vals + par[3] * env_vals^2))
      attract[is.na(attract)] <- 0
      attract <- attract / sum(attract)

      # Sample and update
      step_idx <- sample(length(attract), 1, prob = attract)
      current_cell <- nbhd_cells[step_idx]
      pos <- rowColFromCell(env_raster, current_cell)
      path[i, ] <- c(pos[1], pos[2], current_cell, attract[step_idx])
  }
  path <- as.data.frame(path)
  names(path) <- c("x", "y", "cell", "att")
  return(path)
}
# (Markov version of jag_path is in scratch file, not used for now)

# Calculate variogram of jaguar path
vgram <- function(path, cut = 10, window = 14, start = 1) {
    date0 <- as.Date(path$t_[start])
    date1 <- date0 + window
    path <- path[which(path$t_ >= date0 & path$t_ <= date1), ]
    var <- sapply(1:cut, function(t) {
        p1 <- path[1:(nrow(path) - t), 1:2]
        p2 <- path[(t + 1):nrow(path), 1:2]

        out <- sqrt((p1[, 1] - p2[, 1])^2 + (p1[, 2] - p2[, 2])^2)
        return(mean(unlist(out)))
    })
    plot(var, type = "l")
    return(var)
}

# Plot landscape r with jaguar path and vgram
plot_path <- function(path, int = obs_interval, vgram = FALSE, 
                      type = 1, new = TRUE, ...) {
    # par(mfrow = c(1, ifelse(vgram, 2, 1)))
    # par(mfrow = c(1, 2))
    path <- path[seq(1, nrow(path), int + 1), ]

    col1 <- rgb(1, 0, 0, .5)
    col2 <- rgb(0, 0, 1, .8)
    
    col1 <- magma(nrow(path))

    # Plotting environmental variables + path
    if (new) terra::plot(env_raster)

    if (type == 2) {
      points(path, col = c(col1, col2)[path$state], pch = 19, cex = 0.5)
      for (i in 1:(length(path$state) - 1)) {
          segments(path$x[i], path$y[i], path$x[i + 1], path$y[i + 1], 
                  col = c(col1, col2)[path$state[i]])
      }
      terra::plot(env02[[1]])
      points(path, col = c(col1, col2)[path$state], pch = 19, cex = 0.5)
      for (i in 1:(length(path$state) - 1)) {
          segments(path$x[i], path$y[i], path$x[i + 1], path$y[i + 1], 
                  col = c(col1, col2)[path$state[i]])
      }
    } else if (type == 1) {
      # points(path, col = col1, pch = 19, cex = 0.5)
      points(path, col = rgb(1, 1, 1, 0.3), pch = 19, cex = 0.5)
      lines(path, col = rgb(1, 1, 1, 0.3), lwd = 2)
    }
    # Plotting variogram
    # if (!vgram) return(NULL)
    # plot(vgram(path, ...), type = "l", xlab = "Time lag", ylab = "Variance")
}

# 3. Output analysis -----------------------------------------------------------

results_table <- function(file_ss, file_pp) {
  r_ss <- readRDS(file_ss)
  r_pp <- readRDS(file_pp)
  
  ncol_ss <- length(unlist(r_ss[[1]]))
  ncol_pp <- length(unlist(r_pp[[1]]))
  out_df <- matrix(nrow = nrow(jag_meta), ncol = ncol_ss + ncol_pp + 2)
  for (i in seq_len(nrow(out_df))) {
    if (all(is.na(r_ss[[i]]))) {
      out_df[i, 1:ncol_ss] <- NA
    } else {
      out_df[i, 1:ncol_ss] <- unlist(r_ss[[i]])
      # aic based on likelihood
      out_df[i, ncol_ss + 1] <- 2 * out_df[i, ncol_ss - 1] + 2 * (ncol_ss - 2)
    }
    if (all(is.na(r_pp[[i]]))) {
      out_df[i, (ncol_ss + 2):(ncol_ss + ncol_pp + 1)] <- NA
    } else {
      out_df[i, (ncol_ss + 2):(ncol_ss + ncol_pp + 1)] <- unlist(r_pp[[i]])
      # aic based on likelihood
      out_df[i, ncol_ss + ncol_pp + 2] <- 2 * out_df[i, ncol_ss + ncol_pp] +
                                           2 * (ncol_pp - 2)
    }
  }
  out <- cbind(jag_meta[, c("ID", "biome")], out_df) %>% as.data.frame()
  names(out) <- c("ID", "biome", 
                      paste0("ss_par", 1:9), "ss_ll", "ss_conv", "ss_aic",
                      paste0("pp_par", 1:9), "pp_ll", "pp_conv", "pp_aic")
  return(out)
}

plot_dispersal_comparison <- function(p_ss, p_pp, max_displacement, id) {
                
      plotpdf(nm = paste0("figs/dispersal_tests/indiv_test_", id, "_", 
                          Sys.Date(), ".pdf"), x = 6, y = 6)
      par(mfrow = c(2, 2))

      # Extract probability distributions
      center <- (2 * max_displacement + 1)^2 / 2 + 1
      p_ss[center] <- NA  # Set center to NA for better color scaling
      p_pp[center] <- NA

      # Scale for visualization
      p_ss_scaled <- p_ss / max(p_ss, na.rm = TRUE)
      p_pp_scaled <- p_pp / max(p_pp, na.rm = TRUE)

      # Create rasters
      rast_ss <- rast(matrix(p_ss_scaled, nrow = 2 * max_displacement + 1, 
                             ncol = 2 * max_displacement + 1))
      rast_pp <- rast(matrix(p_pp_scaled, nrow = 2 * max_displacement + 1,
                             ncol = 2 * max_displacement + 1))

      # Define plotting extent based on max_displacement
      init_coords <- xyFromCell(brazil_ras, init_point)
      pixel_res <- res(brazil_ras)[1]
      half_pix <- max_displacement * pixel_res
      plot_extent <- ext(c(
        init_coords[1] - half_pix, init_coords[1] + half_pix,
        init_coords[2] - half_pix, init_coords[2] + half_pix
      ))

      crs(rast_ss) <- crs(brazil_ras)
      ext(rast_ss) <- plot_extent
      crs(rast_pp) <- crs(brazil_ras)
      ext(rast_pp) <- plot_extent

      # Plotting
      track_coords <- self$track[, c("longitude", "latitude")]
      track_vect <- vect(track_coords, geom = c("longitude", "latitude"), crs = wgs84)
      track_vect <- project(track_vect, crs(brazil_ras))
      rast_ss <- terra::crop(rast_ss, track_vect)
      rast_pp <- terra::crop(rast_pp, track_vect)
      terra::plot(rast_ss, main = "Step selection")
      terra::plot(rast_pp, main = "Path propagation")
      terra::plot(rast_ss)
      points(track_vect, pch = 19, cex = 0.3, col = rgb(1, 0, 0, 0.05))
      terra::plot(rast_pp)
      points(track_vect, pch = 19, cex = 0.3, col = rgb(1, 0, 0, 0.05))

      # # Difference plot
      # diff <- p_pp_scaled - p_ss_scaled
      # rast_diff <- rast(matrix(diff, nrow = 2 * max_displacement + 1, ncol = 2 * max_displacement + 1))
      # terra::plot(rast_diff, main = "Difference (PP - SS)")
      # # Forest cover context
      # forest_cover <- terra::crop(brazil_ras[[4]], plot_extent)
      # terra::plot(forest_cover, main = "Forest cover")
      dev.off()
      message("Dispersal comparison plot saved")  
}