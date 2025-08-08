rm(list = ls())
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
# library(plotly)

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

# Brazil biomes shapefile
biome <- vect("data/input/Brazil_biomes/Brazil_biomes.shp")

# Functions ====================================================================

# 0. Basic ---------------------------------------------------------------------

# Logistic function 
exp01 <- function(x) {
    return(1 / (1 + exp(x)))
}

# Calculate moving window average
moving_window <- function(x, y, window = 0.2) {
    out <- sapply(seq_len(length(x)), function(i) {
      pos <- which(x > x[i] - window & x < x[i] + window)
      return(mean(y[pos], na.rm = T))
    })
    df <- data.frame(x = x, y = out)
    return(df[order(df$x), ])
}

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

# Print list of all simulations from simlist.md
simlist <- function() {
    # List all files in the simulations directory
    sim <- readLines("simulations/simlist.md")
    sim <- sim[-which(sim == "")]
    print(sim)
}

# Set up parallel processing
parallel_setup <- function(n_cores = 4) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
    message(paste0("Parallel processing set up with ", n_cores, " cores."))
}

# Load files from a list if they exist in a directory dir
load_if_exists <- function(files, dir) {
    out <- lapply(files, function(f) {
        if (file.exists(paste0(dir, "/", f))) {
            readRDS(paste0(dir, "/", f))
        } else {
            return(NA)
        }
    })
}

# Save raster r as filename fn under ./data/
save_raster <- function(r, fn) {
    # Save raster often because it can get corrupted while working
    writeRaster(r, paste0("data/", fn), overwrite = TRUE)
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


# 1. Data exploration ----------------------------------------------------------

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

# 2. Movement model ============================================================

# Output neighborhood as matrix
# r:    Raster object
# i:    Index of cell in raster
# sz:   Size of neighborhood
make_nbhd <- function(r = brazil_ras, rdf = brdf, i, sz) {
  mat <- matrix(0, nrow = length(i), ncol = (2 * sz + 1)^2)
  # values to add to central cell's row/col to get neighborhood cells' row/col
  ind1 <- t(rep(-sz:sz, each = 2 * sz + 1))
  ind2 <- t(rep(-sz:sz, 2 * sz + 1))
  for (j in seq_len(length(ind1))) {
    mat[, j] <- terra::cellFromRowCol(r, rdf$row[i] + ind1[j], rdf$col[i] + ind2[j])
  }
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
env_function <- function(env, par, nbhd) {
  # Second order
  # attract <- 1 / (1 + exp(par[1] + par[2] * env[, 1] + par[3] * env[, 1]^2 +   # footprint
  #                         par[4] * env[, 2] + par[5] * env[, 2]^2 +            # elevation
  #                         par[6] * env[, 3] + par[7] * env[, 3]^2 +            # slope
  #                         par[8] * env[, 4] + par[9] * env[, 4]^2 +            # forest cover
  #                         par[10] * env[, 5] + par[11] * env[, 5]^2 +          # distance to water
  #                         par[12] * env[, 6] + par[13] * env[, 6]^2))          # distance to road
  
  # # First order no intercept -------------------------------------------------
  # attract <- 1 / (1 + exp(par[1] * env[, 1] + par[2] * env[, 2] +
  #                         par[3] * env[, 3] + par[4] * env[, 4] + 
  #                         par[5] * env[, 5] + par[6] * env[, 6]))
  
  # Simulation -----------------------------------------------------------------
  attract <- 1 / (1 + exp(par[1] + par[2] * env + par[3] * env^2))

  #-----------------------------------------------------------------------------
  attract <- matrix(attract[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  attract <- attract / rowSums(attract, na.rm = TRUE)
  return(attract)  
}

apply_kernel <- function(attract0, kernel) {
  kernel <- kernel / sum(kernel, na.rm = T)
  na_mask <- is.na(attract0)
  attract0[na_mask] <- 0
  p <- attract0 * rep(kernel, each = nrow(attract0))
  p <- p / rowSums(p, na.rm = TRUE)
  p[na_mask] <- NA
  return(p)
}

log_likelihood0 <- function(par, objects, debug = FALSE) {
  nbhd_0   <- objects$nbhd_0
  obs      <- objects$obs
  env_0    <- objects$env_0
  max_dist <- objects$max_dist
  mk       <- objects$mk # only if empirical kernel
  outliers <- objects$outliers

  # Fitted movement kernel -----------------------------------------------------
  # square kernel
  # kernel0 <- dexp(0:(max_dist), rate = exp(par[length(par)])) # move par
  # kernel <- matrix(0, nrow = max_dist * 2 + 1, ncol = max_dist * 2 + 1)
  # center <- max_dist + 1
  # for (i in seq_len(center + max_dist)) {
  #   for (j in seq_len(center + max_dist)) {
  #     kernel[i, j] <- kernel0[max(c(abs(i - center), abs(j - center))) + 1]
  #   }
  # }
  # kernel <- kernel / sum(kernel)
  
  # circular kernel
  # k_exp   <- par[length(par) - 1]
  # bg_rate <- exp01(par[length(par)]) 
  # kernel <- calculate_dispersal_kernel(max_dispersal_dist = max_dist, 
  #                                      kfun = function(x) dexp(x, k_exp) + bg_rate)

  # Empirical movement kernel --------------------------------------------------
  attract0 <- env_function(env_0, par, nbhd_0)
  # attract <- apply_kernel(attract0, kernel)
  attract <- attract0 # no movement kernel
  
  indices <- if (length(outliers) > 0) {
    setdiff(seq_along(obs), outliers)
  } else {
    seq_along(obs)
  }
  like <- sapply(indices, function(i) {
    return(attract[i, obs[i]])
  })

  out <- -sum(log(like), na.rm = TRUE)

  if (is.infinite(out) || is.na(out)) out <- 0
  if (debug) {
    return(list(attract = attract, like = like, out = out, par = par))
  } else {
    return(out)
  }
}

# env_i       : Env variables for each cell of the extended neighborhood
# nbhd_i      : Extended neighborhood for each cell. Rows are path cells, 
#               columns are neighbors.
# to_dest    : Immediate neighbors for each cell. Rows are path cells, 
#              columns are neighbors.
# dest       : 
# obs        : Index of the cell of the extended neighborhood that corresponds
#              to the next GPS observation
# n_obs      : Number of GPS observations
# max_dist   : Maximum distance in pixels for one step
# sim_steps  : Number of steps to simulate
log_likelihood <- function(par, objects, debug = FALSE) {

  env_i      <- objects$env_i
  nbhd_i     <- objects$nbhd_i
  to_dest    <- objects$to_dest
  dest       <- objects$dest
  obs        <- objects$obs
  max_dist   <- objects$max_dist
  sim_steps  <- objects$sim_steps
  outliers   <- objects$outliers
  n_obs      <- length(obs) + 1

  # Movement kernel ------------------------------------------------------------

  # stay/move kernel 
  # move_prob <- exp01(par[length(par)])
  # kernel <- matrix(move_prob/8, nrow = step_size * 2 + 1, ncol = step_size * 2 + 1)
  # center <- step_size + 1
  # kernel[center, center] <- 1 - move_prob

  # square kernel
  # kernel0 <- dexp(0:step_size, rate = par[length(par)]) # move par
  # kernel <- matrix(0, nrow = step_size * 2 + 1, ncol = step_size * 2 + 1)
  # center <- step_size + 1
  # for (i in seq_len(center + step_size)) {
  #   for (j in seq_len(center + step_size)) {
  #     kernel[i, j] <- kernel0[max(c(abs(i - center), abs(j - center))) + 1]
  #   }
  # }
  # kernel <- kernel / sum(kernel)

  # circular kernel
  # k_exp   <- par[length(par) - 1] # second-to-last parameter
  # bg_rate <- exp01(par[length(par)]) # background rate
  # kernel <- calculate_dispersal_kernel(max_dispersal_dist = step_size, 
  #                                      kfun = function(x) dexp(x, k_exp) + bg_rate)

  # Attractiveness function ----------------------------------------------------
  attract0 <- env_function(env_i, par, nbhd = nbhd_i) 
  attract  <- attract0   # no movement kernel
  # attract <- apply_kernel(attract0, kernel)

  # Propagating probabilities forward ------------------------------------------
  ncell_local <- (2 * max_dist + 1)^2 
  current <- array(0, dim = c(ncell_local, n_obs, sim_steps))
  
  # center     : Center of ncell_local (center cell of (2 * buffer + 1)
  # Set to 1 at step #1 for each observation because that's where it actually is
  center <- ncell_local / 2 + 0.5
  current[center, , 1] <- 1

  for (j in 1:(sim_steps - 1)) {
    # Probabilities across entire nbhd_i for step j
    step_prob <- as.vector(current[, , j]) * attract[]
    #   e.g. if to_dest[1, 2] = 3, then the 2nd neighbor of the 1st cell of
    #   nbhd_i is the 3rd cell of nbhd_i.
    dest[] <- step_prob[as.vector(to_dest)]
    # Summing probabilities up to step j to generate step j+1
    current[, , j + 1] <- rowSums(dest, na.rm = TRUE)
  }

  # Calculate log likelihood 
  predictions <- matrix(0, nrow = sim_steps, ncol = n_obs)
  for (i in 1:n_obs) {
    if (i %in% outliers) {
      predictions[, i] <- rep(NA, sim_steps)
      next
    }
    prob <- current[obs[i], i, ]
    predictions[, i] <- prob #+ bg_rate - prob * bg_rate
  }

  log_likelihood <- rowSums(log(predictions), na.rm = TRUE) 
  # out            <- -max(log_likelihood, na.rm = TRUE)
  out <- -log_likelihood[sim_steps]
  
  if (is.infinite(out) || is.na(out)) out <- 0
  if (debug) {
    return(list(out = out, array = current, predictions = predictions, par = par))
  } else {
    return(out)
  }
}

# 2b. Simulation ---------------------------------------------------------------

# Generates a random field with a given correlation structure
# b: beta parameter for gstat, s: sill, r: range, n: nugget
gen_landscape <- function(size = 100, b = 1, s = 0.03, r = 10, n = 0) {   
    xy <- expand.grid(1:size, 1:size)
    names(xy) <- c("x", "y")

    # Autocorrelation model
    model <- gstat(formula = z ~ 1, locations = ~x + y, dummy = TRUE, beta = b, 
                   model = vgm(psill = s, range = r, model = "Exp", nugget = n), 
                   nmax = 20)
    out <- predict(model, newdata = xy, nsim = 1)
    if (any(out < 0)) out[out < 0] <- 0

    # Output as raster 
    out <- rast(out)
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

  # Movement kernel (circular)
  # k_exp   <- par[length(par) - 1]
  # bg_rate <- exp01(par[length(par)]) 
  # kern_m <- calculate_dispersal_kernel(max_dispersal_dist = step_size, 
  #                                      kfun = function(x) dexp(x, k_exp) + bg_rate)

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

## Load parameters and likelihood
load_output <- function(name, i = NULL) {
    dir <- paste0("data/output/sim_", name, "/")
    # ll_files <- list.files(dir)[grep("likelihood_", list.files(dir))]
    # par_files <- list.files(dir)[grep("par_out_", list.files(dir))]
    # ll_files <- paste0("likelihood_", 1:njag, ".rds")
    # par_files <- paste0("par_out_", 1:njag, ".rds")
    seq <- if (is.null(i)) 1:njag else i
    out_files <- paste0("out_", seq, ".rds")
    out <- load_if_exists(out_files, dir)
    if (length(out) == 1) out <- out[[1]]
    return(out)
}

par_to_df <- function(par) {
    df <- do.call(rbind, lapply(par, function(x) {
        print(x[[1]])
    }))
}

results_table <- function(s, params = TRUE, holdout = TRUE) {

 npar <- vector()
  
 lldf <- sapply(s, function(i) {
    indiv <- paste0("out_", 1:82, ".rds")
    out <- sapply(indiv, function(j) {
      print(paste(i, j))
      if (!(j %in% list.files(paste0("data/output/sim_", i, "/")))) {
        return(NA)
      } else {
        out <- readRDS(paste0("data/output/sim_", i, "/", j))
        npar[i] <<- length(out$par)
        # return(c(out$out, out$par))
        return(out$out)
      }
    })
    return(out)
  }) %>% data.table()
  colnames(lldf) <- paste0("ll_", s)

  aic <- lldf * 2 + npar[col(lldf)] * 2 
  colnames(aic) <- paste0("aic_", s)
  df <- cbind(lldf, aic)

  df <- cbind(jag_id, sex = as.factor(jag_meta$Sex),
                  age = as.numeric(jag_meta$Estimated.Age),
                  weight = as.numeric(jag_meta$Weight),
                  bio = as.factor(jag_meta$biome), df,
                  nmove = as.numeric(jag_meta$nmove),
                  ndays = as.numeric(jag_meta$ndays),
                  meandist = as.numeric(jag_meta$meandist),
                  totdist = as.numeric(jag_meta$totdist))

  return(df)
}

# Given fitted parameter values and individual id, plot local environment transformed
# by environmental function.
plot_attraction <- function(id, results, model = 2) {
  id <- as.numeric(id)
  pars <- switch(model, results[results$id == id, 1:6], results[results$id == id, 9:14])
  # make neighborhood for individual (extended)
  path <- jag_move[ID == id]
  # get bounding box of path and crop brazil_ras
  gps <- vect(path[, 3:4], geom = c("longitude", "latitude"))
  # add small buffer to bbox
  bbox <- ext(gps) #+ c(-0.1, 0.1, -0.1, 0.1) # 0.1 degrees buffer
  # get raster of local environment
  local_ras <- crop(brazil_ras, bbox)
  # normalize each layer of local_ras to [-1, 1]
  local_ras <- scale(local_ras, center = TRUE, scale = TRUE)

  attract <- 1 / (1 + exp(pars[[1]] * local_ras[[1]] +
                     pars[[2]] * local_ras[[2]] + 
                     pars[[3]] * local_ras[[3]] + 
                     pars[[4]] * local_ras[[4]] + 
                     pars[[5]] * local_ras[[5]] + 
                     pars[[6]] * local_ras[[6]]))
                     
  # plot local environment in terms of attraction
  terra::plot(attract, main = paste("Local environment for jaguar", id))
  # plot jaguar path points as arrows between points
  
  n <- nrow(path)
  
  colors <- colorRampPalette(c("gold", "darkorange", "darkred"))(n-1)

  points(path$longitude, path$latitude, pch = 19, cex = 0.5, col = rgb(1, 0, 0, 0.5))
  # arrows(x0 = path$longitude[-n],
  #       y0 = path$latitude[-n],
  #       x1 = path$longitude[-1], 
  #       y1 = path$latitude[-1],
  #       length = 0.08,
  #       col = colors,
  #       lwd = 1.5)
}
  # plot local environment in terms of attraction
  # terra::plot(attract, main = paste("Local environment for jaguar", id))
  # plot jaguar path points as arrows between points
  
  # n <- nrow(path)
  # colors <- colorRampPalette(c("gold", "darkorange", "darkred"))(n-1)
  # arrows(x0 = path$longitude[-n],
  #       y0 = path$latitude[-n],
  #       x1 = path$longitude[-1], 
  #       y1 = path$latitude[-1],
  #       length = 0.08,
  #       col = colors,
  #       lwd = 1.5)