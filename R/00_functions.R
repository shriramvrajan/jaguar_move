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
library(fitdistrplus)
library(foreach)
library(doParallel)
library(viridis)
# library(plotly)

# Global parameters ============================================================

true_step <- 1     # How many pixels does jaguar move in 1 time step?

# CRS definitions
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
epsg5880 <- "+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 
+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"

# Simulation parameters
# sim_n <- 1000

# Data =========================================================================

# Jaguar movement data, ID numbers, and metadata
jag_move <- readRDS("data/jag_data_BR.rds")
jag_id <- readRDS("data/jag_list.rds")
njag <- nrow(jag_id)
jag_meta <- data.table(read.csv("data/input/jaguars/jaguar_metadata2.csv"))

# RasterStack of environmental variables -- see 01_generate_data.R for details
brazil_ras <- rast("data/env_layers.grd")
load("data/env_layers.RData") # same information as a dataframe ('brdf')

# Brazil biomes shapefile
biome <- vect("data/input/Brazil_biomes/Brazil_biomes.shp")

# Functions ====================================================================

# 0. Basic ---------------------------------------------------------------------

# Logistic function 
exp01 <- function(x) {
    return(exp(x) / (1 + exp(x)))
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

# Set up parallel processing
parallel_setup <- function(n_cores = 4) {
    cl <- makeCluster(n_cores)
    registerDoParallel(cl)
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
    names(outdf)[2:3] <- c("row", "col")
    return(outdf)
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

# Produce amt::track object for jaguar i
make_track0 <- function(id) {
    id <- as.numeric(id)
    path <- jag_move[ID == id]
    path$t <- lubridate::mdy_hm(as.character(path$timestamp))
    path <- vect(path, geom = c("longitude", "latitude"), crs = wgs84)
    path <- project(path, epsg5880)
    path <- track(x = crds(path)[, 1], y = crds(path)[, 2], 
                  t = path$t, id = path$ID, crs = epsg5880)
}

# Decompose timestamp into year, month, day, hour and add metadata to track
make_full_track <- function(id) {
    dat <- jag_move[ID == id]
    dat$year <- as.numeric(format(as.POSIXct(dat$timestamp, 
                            format = "%m/%d/%Y %H:%M"), "%Y"))
    dat$year <- ifelse(dat$year > 23, dat$year + 1900, dat$year + 2000)

    dat$mon <- as.numeric(format(as.POSIXct(dat$timestamp, 
                           format = "%m/%d/%Y %H:%M"), "%m"))
    dat$day <- as.numeric(format(as.POSIXct(dat$timestamp, 
                           format = "%m/%d/%Y %H:%M"), "%d"))
    dat$hr <- format(as.POSIXct(dat$timestamp, 
                           format = "%m/%d/%Y %H:%M"), "%H:%M")
    dat$hr <- as.numeric(gsub(":[0-9][0-9]", "", dat$hr))
    dat <- dat[, timestamp := NULL]
    tr <- make_track0(id)
    st <- steps(tr)
    dat$sl <- c(NA, st$sl_)             # step lengths in m
    dat$ta <- c(NA, st$ta_)             # turn angles in radians
    dat$dir <- c(NA, st$direction_p)    # bearing in radians
    dat$dt <- c(NA, as.numeric(st$dt_)) # time interval in minutes
    dat$spd <- dat$sl / dat$dt
    return(dat[, c("longitude", "latitude", "ID", "year", "mon", "day", "hr", 
                   "sl", "ta", "dir", "dt", "spd")])
}

# 1. Data exploration ----------------------------------------------------------

# Plot environmental layers for path 
# layer: layer number of raster
# bounds: bounding box of path
# path: path of jaguar
# grad: number of colors to use for path
# ras: raster stack of environmental layers
plot_env <- function(layer, bounds, path, grad = 20, ras = env) {
    terra::plot(crop(ras[[layer]], bounds), main = names(ras)[layer], 
                col = (gray.colors(50, start = 0.7, end = 0.9)))
    # pal <- colorRampPalette(c("#ff0000", "#1bbe29"))(grad)
    pal <- rep(rgb(0.5, 0.5, 0.5, 0.5), 20)
    path <- as.data.frame(path)
    path$col <- rep(pal, each = ceiling(nrow(path) / grad))[seq_len(nrow(path))]
    points(x = path$x, y = path$y, pch = 19, cex = 0.8, col = path$col)
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
map_track <- function(i, grad = 20, type = 2) {
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
plot_the_curve <- function(par, bounds = c(0, 10)) {
    # Plot functional form for movement model
    x <- seq(bounds[1], bounds[2], 0.1)
    y <- 1 / (1 + exp(par[1] + par[2] * x + par[3] * x^2))
    plot(x, y, type = "l", col = "blue", lwd = 3, xlab = "Environmental variable", 
         ylab = "Attractiveness")
    abline(v = 0, lty = 2)
    abline(h = 0.5, lty = 2)
}

### 2. Movement model ==========================================================

# Output neighborhood as matrix ------------------------------------------------
# r:    Raster object
# rdf:  Dataframe of raster cells
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

# Normalize probabilities across neighbors of each cell ------------------------
# v: vector of cell values
normalize_nbhd <- function(v) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = TRUE)
}

# Prepare input objects for movement model -------------------------------------
# traject:   Path of individual, vector of raster cell indices
# max_dist:  Maximum distance in pixels for one step
# nbhd0:     Global neighborhood 
# r:         Raster object
# rdf:       Dataframe of raster cells
prep_model_objects <- function(traject, max_dist, r) {
    
    # Make global neighborhood from raster
    message("Building global neighborhood")
    rdf <- raster_to_df(r)
    nbhd0 <<- make_nbhd(i = seq_len(nrow(rdf)), sz = true_step, 
                        r = r, rdf = rdf) 

    # Extended neighborhoods of each cell in individual's trajectory
    message("Building neighborhoods for each cell")
    nbhd_index <- make_nbhd(i = traject, sz = max_dist, r = r, rdf = rdf)

    # Each entry in the list is the immediate neighborhood of each cell in the 
    # extended neighborhood, as represented by a cell number of raster r
    message("Getting indices of extended neighborhood of each cell")
    nbhd_list <- lapply(seq_len(nrow(nbhd_index)), function(i) {                
      row_inds <- seq_len(ncol(nbhd_index)) + (i - 1) * ncol(nbhd_index)
      names(row_inds) <- nbhd_index[i, ] # cell numbers as names for indexing
      out <- matrix(row_inds[as.character(nbhd0[nbhd_index[i, ], ])], 
                    nrow = length(row_inds), ncol = ncol(nbhd0))
      return(out)
    })
    nbhd <<- do.call(rbind, nbhd_list)
    # Reindexing allows linkage of row numbers from nbhd to raster cells
    nbhd_index <- as.vector(t(nbhd_index))
    nbhd_index <<- nbhd_index
    
    message("Getting indices of immediate neighborhood of each cell")
    # For each cell of the extended neighborhood of the path, what are
    # the immediate neighbors? Rows are each cell of nbhd, columns are row #s 
    # from nbhd. All row lengths standardized with missing neighbors as NAs.
    to_dest <- tapply(seq_len(length(nbhd)), nbhd, function(x) {  
      if (length(x) > 9) print(x)
      out <- c(x, rep(NA, ncol(nbhd) - length(x)))
      return(out)
    })
    to_dest <<- t(matrix(unlist(to_dest), nrow = ncol(nbhd), ncol = nrow(nbhd)))
    dest <<- matrix(0, nrow = nrow(nbhd), ncol = ncol(nbhd))

    message("Indexing observed data")
    # Building observed data to test against
    index_mat <- matrix(
      data = seq_len(length(nbhd_index)),
      nrow = (nrow(nbhd) / length(traject)),
      ncol = length(traject)
    )
    index_mat <<- index_mat
    # For each step, which cell of the extended nbhd did it go to next?
    obs <- vector(length = ncol(index_mat) - 1)
    for (y in 1:(ncol(index_mat) - 1)) {                                         # 13s
      test <- which(nbhd_index == traject[y + 1])
      num <- which(index_mat[1, y] <= test & test <= index_mat[nrow(index_mat), y])
      obs[y] <- which(index_mat[, y] == test[num])
    } 
    obs <<- obs

    message("Getting environmental variables")
    # Normalizing desired environmental variables for extended neighborhood    
    # env <- scales::rescale(rdf$sim1[nbhd_index], to = c(0, 1))  
    env <- rdf$sim1[nbhd_index]
    names(env) <- seq_len(length(nbhd_index))
    env <<- env

    message("Prepared model objects.")
}

# Attractiveness function for movement model -----------------------------------
# env:      Environmental variable values
# par:      Parameters for functional form
# format:   TRUE to return as matrix, FALSE to return as vector
env_function <- function(env, par, format = FALSE) {
  attract <- 1 / (1 + exp(par[1] + par[2] * env + par[3] * env^2))
  if (format) {
    attract <- matrix(attract[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  }
  return(attract)
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

make_movement_kernel1 <- function(n = 10000, sl_emp, ta_emp, max_dist, minimum = 0) {
  
  # Fit gamma parameters? 
  avail <- data.frame(sl = sample(sl_emp, n, replace = TRUE),
                      ta = sample(ta_emp, n, replace = TRUE))
  avail$x <- avail$sl * cos(avail$ta) / 1000
  avail$y <- avail$sl * sin(avail$ta) / 1000
  avail$d <- floor(sqrt(avail$x^2 + avail$y^2))

  stay_prob <- length(which(avail$d == 0)) / length(avail$d)  

  moves <- avail$d[-which(avail$d == 0)]
  rate <- fitdist(moves, distr = "exp", method = "mle")$estimate
  # Return out$n as a vector of probabilities
  # But using 0-1 for the center cell + gamma dist for the rest

  # OK YOU CAN DO THIS GANBATTE !!!!!!!!!!!!!

  out <- data.frame(x = rep(-max_dist:max_dist, times = size),
                    y = rep(-max_dist:max_dist, each = size))
  return(out$n)

  # avail$xi <- sapply(avail$x, function(x) ifelse(x > 0, floor(x), ceiling(x)))
  # avail$yi <- sapply(avail$y, function(y) ifelse(y > 0, floor(y), ceiling(y)))
  # density <- tapply(avail$x, list(avail$yi, avail$xi), length)
  # density <- reshape2::melt(density)
  # names(density) <- c("x", "y", "n")
  # size <- max_dist * 2 + 1

  # out <- merge(out, density, by = c("x", "y"), all.x = TRUE)
  # out$n[is.na(out$n)] <- 0
  # out$n <- out$n + minimum
  # out$n <- out$n / sum(out$n)  
}



# Return -(maximum log likelihood) given a set of parameters -------------------
# For traditional SSF (log_likelihood0) and all others (log_likelihood)
# par:      Parameters for optimization
# objects:  List of objects for optimization
log_likelihood0 <- function(par, objects) {
  # par        : Initial values of parameters for optimization
  env <- objects[[1]]
  max_dist <- objects[[2]]
  mk <- objects[[3]]
  obs <- objects[[4]]
  n_obs <- length(obs) + 1
  
  # Attractiveness function 0: traditional SSF 
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])

  # Attractiveness function 1: simulations
  attract_e <- normalize_nbhd(env_function(env, par[2:4]))

  step_range <- (max_dist * 2 + 1) ^ 2
  
  p_obs <- sapply(seq_len(n_obs - 1), function(t) {
    env_local <- attract_e[(step_range * (t - 1) + 1):(step_range * t)]
    env_local <- env_local / sum(env_local, na.rm = TRUE)
    # TESTING WITHOUT MK -------------------------------------------------------
    # p <- env_local
    # With MK ------------------------------------------------------------------
    env_weight <- exp01(par[1])
    p <- env_weight * env_local + (1 - env_weight) * mk # mk = movement kernel
    return(ifelse((p[obs[t]] == 0), 0.001, p[obs[t]]))
  })
  ll <- -sum(log(p_obs))
  if (is.infinite(ll)) ll <- 0
  # print(ll)
  return(ll)
}

log_likelihood <- function(par, objects, debug1 = FALSE) {
  # Environmental variables
  env       <- objects[[1]]
  # Neighborhood
  nbhd       <- objects[[2]]
  # Maximum distance in pixels for one step
  max_dist   <- objects[[3]]
  # Number of GPS observations (length of track)
  sim_steps  <- objects[[4]]
  # to_dest    : For each cell of the extended neighborhood of the path, what 
  #               are the immediate neighbors? Rows are path cells, columns are 
  #               neighbors.
  to_dest    <- objects[[5]]
  # obs        : Index of the cell of the extended neighborhood that corresponds
  #              to the next GPS observation
  obs        <- objects[[6]]
  n_obs      <- length(obs) + 1

  # Attraction function 3: simulations -----------------------------------------
  # attract1 <- normalize_nbhd(exp(par[1] * env)) # + exp(par[2] * env2)
  # print(par)
  # attract1 <- env_function(env, par[2:4])
  # move_prob <- exp01(par[1])
  # attract <- t(apply(attract1, 1, function(r) {
  #   cent <- ceiling(length(r) / 2)
  #   r[cent] <- r[cent] * (1 - move_prob)
  #   r[-cent] <- r[-cent] * ((move_prob) / (sum(!is.na(r)) - 1))
  #   return(r / sum(r, na.rm = TRUE))
  # }))

  # Attraction function 3a: simulation, no move param --------------------------
  attract <- normalize_nbhd(env_function(env, par))

  # Attraction function 4: With 0-1 parameter ----------------------------------
  # move_prob <- exp01(par[7])
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  # # attract_h <- exp(par[8] * env$home)
  # # attract <- normalize_nbhd(attract_e * attract_h)
  # attract <- normalize_nbhd(attract_e)
  # attract <- t(apply(attract, 1, function(r) {
  #   cent <- ceiling(length(r) / 2)
  #   r[cent] <- r[cent] * (1 - move_prob)
  #   r[-cent] <- r[-cent] * (move_prob / (sum(!is.na(r)) - 1))
  #   return(r / sum(r, na.rm = TRUE))
  # }))

  # Array for propagating probabilities forward ================================
  # n_obs      : Number of GPS observations
  # steps      : Number of simulated steps
  step_range <- (2 * max_dist + 1)^2 
  current <- array(0, dim = c(step_range, n_obs, sim_steps))
  
  # center     : Center of step_range (center cell of (2 * buffer + 1)
  # Set to 1 at step #1 for each observation because that's where it actually is
  center <- step_range / 2 + 0.5
  current[center, , 1] <- 1

  for (j in 1:(sim_steps - 1)) {
    # Probabilities across entire nbhd for step j
    step_prob <- as.vector(current[, , j]) * attract[]

    # dest has same dimensions as nbhd
    # to_dest is a matrix with the same dimensions as nbhd
    #   rows are neighborhood cells, columns are neighbors of those cells
    #   each entry is the index of the neighbor cell
    #   e.g. if to_dest[1, 2] = 3, then the 2nd neighbor of the 1st cell of
    #   nbhd is the 3rd cell of nbhd.
    dest[] <- step_prob[as.vector(to_dest)]
    # Summing probabilities up to step j to generate step j+1
    current[, , j + 1] <- rowSums(dest, na.rm = TRUE)
  }
  
  #### Run likelihood function from fitting with true parameters

  # Calculate log likelihood 
  predictions <- matrix(0, nrow = sim_steps, ncol = n_obs)
  for (i in 1:n_obs) {
    prob <- current[obs[i], i, ]
    predictions[, i] <- ifelse(prob == 0, 0.001, prob)
    # returns the probability for the row associated with the next 
    # observation location, for that observation i, across all time stepsp
  }

  log_likelihood <- rowSums(log(predictions), na.rm = TRUE)
  # log of product is sum of logs

  # out <- -max(rowSums(log(predictions), na.rm = TRUE), na.rm = TRUE)
  ## DEBUG
  out <- -log_likelihood[sim_interval + 2]
  if (is.infinite(out) || is.na(out)) out <- 0

  if (debug1) {
    return(list(out = out, predictions = predictions))
  } else {
    return(out)
  }
  # Return negative of the maximum log likelihood because we want to minimize
  # Lower negative log likelihood = higher likelihood 
}

loglike <- function(par, objects) {
    # Wrapper function for log_likelihood and log_likelihood0
    if (refit_model0 == TRUE) {
        # Traditional SSF: 
        # Intermediate spaces between GPS observations unimportant
        return(log_likelihood0(par, objects))
    } else {
        # New model: intermediate spaces important
        return(log_likelihood(par, objects))
    }
}

run_optim <- function(param, objects, i) {
    ntries <- 0
      ## Main fitting loop, tries each individual 20x and moves on if no fit
    while (ntries <= 20) {
        tryCatch({
            par_out <- optim(param, loglike, objects = objects)
            saveRDS(par_out, paste0("par_out_", i, ".rds"))

            message("Running loglike_fun...")
            likelihood <- loglike(par_out[[1]], objects = objects)
            saveRDS(likelihood, paste0("likelihood_", i, ".rds"))

            message(paste0("jaguar ", i, " fitted ", date()))
            ntries <- 21 # End while loop
          },
          error = function(e) {
            message(e)
            message(paste("Try #:", ntries))
            if (ntries == 20) {
              message("Skipping, couldn't fit in 20 tries")
              saveRDS(NA, paste0("NA_", i, ".rds"))
            } else {
              message("Retrying")
            }
          },
          finally = {
            ntries <- ntries + 1
          }
        )
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

    # Output as both raster and data frame
    out <- rast(out)
    # plot(out)
    outdf <- raster_to_df(out)
    return(list(raster = out, df = outdf))
}

# Generate a jaguar path of n steps starting from (x0, y0) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x0, y0, n_step, par, neighb) {
    # Set initial state
    path <- matrix(NA, nrow = n_step, ncol = 4)
    path[1, ] <- c(x0, y0, NA, NA)

    if (length(par) == 3) par <- c(0, par)

    # Probability of moving from current grid cell
    move_prob <- exp01(par[1])

    for (i in 2:n_step) {
        if (i %% 10 == 0) print(i)
        pos <- path[i - 1, 1:2]

        # Attractiveness of each cell in neighborhood
        nbhd <- as.vector(make_nbhd(r = env01[[1]], rdf = env01[[2]], sz = neighb,
                          i = terra::cellFromRowCol(env01[[1]], pos[1], pos[2])))
        att <- sapply(nbhd, function(x) {
          if (is.na(x)) {
            return(NA)
          } else {
            return(env_function(env01[[1]][x], par[2:4], format = FALSE)$sim1)
          }
        })
        if (any(is.na(att))) att[is.na(att)] <- 0
        att <- att / sum(att)
        attract <- att

        # Find center cell and scale by move_prob (COMMENT OUT FOR NO MOVEPAR)
        # cent <- ceiling(length(att) / 2)
        # att[cent] <- att[cent] * (1 - move_prob)
        # att[-cent] <- att[-cent] * (move_prob / (sum(!is.na(att)) - 1))
        # attract <- att / sum(att)

        # Sample next step
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- c(rowColFromCell(env01[[1]], nbhd[step]), 
                       nbhd[step], attract[step])
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
plot_path <- function(path, int = sim_interval, vgram = FALSE, 
                      type = 1, new = TRUE, ...) {
    # par(mfrow = c(1, ifelse(vgram, 2, 1)))
    # par(mfrow = c(1, 2))
    path <- path[seq(1, nrow(path), int + 1), ]

    col1 <- rgb(1, 0, 0, .5)
    col2 <- rgb(0, 0, 1, .8)
    
    col1 <- magma(nrow(path))

    # Plotting environmental variables + path
    if (new) terra::plot(env01[[1]])

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
      points(path, col = "black", pch = 19, cex = 1)
      lines(path, col = rgb(0, 0, 0, 0.5), lwd = 2)
    }
    # Plotting variogram
    # if (!vgram) return(NULL)
    # plot(vgram(path, ...), type = "l", xlab = "Time lag", ylab = "Variance")
}

# 3. Output analysis -----------------------------------------------------------

## Load parameters and likelihood
load_output <- function(name) {
    dir <- paste0("data/output/", name, "/")
    # ll_files <- list.files(dir)[grep("likelihood_", list.files(dir))]
    # par_files <- list.files(dir)[grep("par_out_", list.files(dir))]
    ll_files <- paste0("likelihood_", 1:njag, ".rds")
    par_files <- paste0("par_out_", 1:njag, ".rds")
    ll <- load_if_exists(ll_files, dir)
    par <- load_if_exists(par_files, dir)
    return(list(unlist(ll), par))
}

par_to_df <- function(par) {
    df <- do.call(rbind, lapply(par, function(x) {
        print(x[[1]])
    }))
}

results_table <- function(s = c("K", "RW", "RWM", "RWH", "trad2")) {
  
  # null <- load_output("LL_null")
  
  output <- lapply(s, function(x) {
    load_output(paste0("LL_", x))
  })

  likelihood <- lapply(output, function(x) -x[[1]])
  params <- lapply(output, function(x) par_to_df(x[[2]]))
  names(params) <- s
  npar <- sapply(params, function(x) ncol(x))
  aic <- lapply(seq_len(length(likelihood)), function(i) {
    2 * npar[i] - 2 * likelihood[[i]]
  })  
  

  df <- as.data.frame(cbind(do.call(cbind, likelihood), do.call(cbind, aic)))
  colnames(df) <- c(paste0("l_", s), paste0("aic_", s))

  ## Extra variables
  nmove <- sapply(1:njag, function(x) {
    length(which(jag_move$ID == as.numeric(jag_id[x])))
  })
  ndays <- sapply(1:njag, function(x) {
    moves <- jag_move[ID == as.numeric(jag_id[x])]
    dates <- sort(as.Date(sapply(moves$timestamp, function(dt) {
            strsplit(as.character(dt), " ")[[1]][1]
        }), format = "%m/%d/%y"))
    return(as.numeric(difftime(dates[length(dates)], dates[1])))
  })
  meandist <- tapply(jag_move$dist, jag_move$ID, function(x) mean(x, na.rm = TRUE))
  totdist <- tapply(jag_move$dist, jag_move$ID, function(x) sum(x, na.rm = TRUE))

  df <- cbind(jag_id, sex = as.factor(jagmeta_br$Sex),
                  age = as.numeric(jagmeta_br$Estimated.Age),
                  weight = as.numeric(jagmeta_br$Weight),
                  bio = as.factor(jagmeta_br$biome),
                  nmove, ndays, meandist, totdist, df)

  return(list(df, params))
}

### EXTRA LOG LIKELIHOOD METHODS ===============================================

  # Attraction function 1: environmental variables + home range ----------------
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6] +
  #                  par[7] * env$home)
  # attract_h <- exp(par[7] * env$home)
  # attract <- normalize_nbhd(attract_e * attract_h) #* normalize_nbhd(attract_t)
  # attract <- normalize_nbhd(attract_e)

  # Attraction function 2: just home range -------------------------------------
  # attract_h <- exp(par[1] * env$home)
  # attract <- normalize_nbhd(attract_h) 