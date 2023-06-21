rm(list = ls())

## Basic functions and libraries

# Libraries
library(terra)
library(tidyverse)
library(data.table)
library(finch)
library(gstat)
library(ctmm)
library(amt) 
library(lubridate)
# library(mixtools)

# Output message m both in console and in logfile f
msg <- function(m, f = "data/output/run_log.txt") {
    m <- paste(format(Sys.time(), "%d.%m.%y  %R"), m)
    print(m)
    cat(m, file = f, append = TRUE, sep = "\n")
}

# Global parameters ============================================================

buffersize <- 1 # How far does jaguar move in 1 time step

wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +type=crs"
epsg5880 <- "+proj=poly +lat_0=0 +lon_0=-54 +x_0=5000000 +y_0=10000000 
+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs"


# Functions ====================================================================

# Basic ------------------------------------------------------------------------

# Save raster x as filename fn under ./data/
save_ras <- function(x, fn) {
    # Save raster often because it can get corrupted while working
    writeRaster(x, paste0("data/", fn), overwrite = TRUE)
}

rast_df <- function(r) {
    outdf <- as.data.frame(r)
    outdf <- cbind(outdf, rowColFromCell(r, seq_len(nrow(outdf))))
    return(outdf)
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

# Data exploration -------------------------------------------------------------

# Plot environmental layers for path 
plot_env <- function(layer, bounds, path, ras = env) {
    terra::plot(crop(ras[[layer]], bounds), main = names(ras)[layer])
    lines(as.data.frame(path), col = rgb(1, 0, 0, 0.3), pch = 4)
}

# Map path of jaguar i
map_jag <- function(i, type = 2) {
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
            plot_env(layer = i, bounds = bbox, path = path)
        }
    }
}

# Map home range of jaguar i
map_homerange <- function(id, vgram = FALSE) {
  msg(paste0("Mapping home range of jaguar ", id, "..."))
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

# Produce amt::track object for jaguar i
jag_track <- function(id) {
    id <- as.numeric(id)
    path <- jag_move[ID == id]
    path$t <- lubridate::mdy_hm(as.character(path$timestamp))
    path <- vect(path, geom = c("longitude", "latitude"), crs = wgs84)
    path <- project(path, epsg5880)
    path <- track(x = crds(path)[, 1], y = crds(path)[, 2], 
                  t = path$t, id = path$ID, crs = epsg5880)
}

# Movement model ---------------------------------------------------------------

# Outputs a matrix of cell numbers corresponding to raster (r, rdf)
# based on a central cell (i) and a buffer size around that cell (sz)
make_nbhd <- function(r = brazil_ras, rdf = brdf, i, sz) {
  # if x is the center square and o are neighbors, and e.g. sz = 2
  # (2*sz + 1)^2 represents total neighborhood size 
  mat <- matrix(0, nrow = length(i), ncol = (2 * sz + 1)^2)
  # values to add to central cell's row/col to get neighborhood cells' row/col
  ind1 <- t(rep(-sz:sz, each = 2 * sz + 1))
  ind2 <- t(rep(-sz:sz, 2 * sz + 1))
  for (j in seq_len(length(ind1))) {
    mat[, j] <- cellFromRowCol(
      r, rdf[i, 2] + ind1[j],
      rdf[i, 3] + ind2[j]
    )
  }
  return(mat)
}

# Normalize probabilities across neighbors of each cell
norm_nbhd <- function(v) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = TRUE)
}

# Prepare input objects for movement model
input_prep <- function(traject, max_dist, steps, nbhd0, r, rdf) {

    # Extended neighborhoods of each cell in individual's trajectory
    msg("Building neighborhoods for each cell")
    nbhd_index <- make_nbhd(i = traject, sz = max_dist, r = r, rdf = rdf)
  
    # Each entry in the list is the immediate neighborhood of each cell in the 
    # extended neighborhood, as represented by a cell number of raster r
    msg("Getting indices of extended neighborhood of each cell")
    nbhd_list <- lapply(seq_len(nrow(nbhd_index)), function(i) {                 # 14s
      # msg(i)
      # For each actual cell in the path, what are the cells in the extended
      # neighborhood? Rows are path cells, columns are extended neighborhood.
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

    msg("Getting indices of immediate neighborhood of each cell")
    # For each cell of the extended neighborhood of the path, what are
    # the immediate neighbors? Rows are path cells, columns are neighbors.
    # All row lengths standardized by turning missing neighbors into NAs.
    to_dest <- tapply(seq_len(length(nbhd)), nbhd, function(x) {                 # 36s
      # print(nbhd[x])
      # print(ncol(nbhd) - length(x))
      out <- c(x, rep(NA, ncol(nbhd) - length(x)))
      return(out)
    })
    to_dest <<- t(matrix(unlist(to_dest), nrow = ncol(nbhd), ncol = nrow(nbhd)))
    dest <<- matrix(0, nrow = nrow(nbhd), ncol = ncol(nbhd))

    msg("Indexing observed data...")
    # Building observed data to test against
    index_mat <- matrix(
      data = seq_len(length(nbhd_index)),
      nrow = (nrow(nbhd) / length(traject)),
      ncol = length(traject)
    )
    index_mat <<- index_mat
    obs <- vector(length = ncol(index_mat) - 1)
    for (y in 1:(ncol(index_mat) - 1)) {                                         # 13s
      test <- which(nbhd_index == traject[y + 1])
      num <- which(index_mat[1, y] < test & test < index_mat[nrow(index_mat), y])
      obs[y] <- which(index_mat[, y] == test[num])
    } 
    obs <<- obs
}

# Returns negative of the maximum log likelihood given a set of parameters (par)  
loglike_fun <- function(par) {
  # par        : Initial values of parameters for optimization
  # nbhd       : Neighborhood
  # step_range : Step range
  # n_obs      : Number of GPS observations (length of track)
  # steps:     : Number of steps simulated
  # to_dest    : For each cell of the extended neighborhood of the path, what 
  #               are the immediate neighbors? Rows are path cells, columns are 
  #               neighbors.
  # obs        : Index of the cell of the extended neighborhood that corresponds
  #               to the next GPS observation
  # env        : Environmental variables

  # Attractiveness function 1: environmental variables + home range ------------
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  # attract_h <- exp(par[7] * env$home)
  # attract <- norm_nbhd(attract_e * attract_h) #* norm_nbhd(attract_t)

  # Attractiveness function 2: just home range ---------------------------------
  # attract_h <- exp(par[1] * env$home)
  # attract <- norm_nbhd(attract_h) 

  # Attractiveness function 3: simulations -------------------------------------
  attract1 <- norm_nbhd(exp(par[1] * env1)) # + exp(par[2] * env2)
  attract2 <- norm_nbhd(exp(par[2] * env2))
  attract <- attract1

  # Array for propagating probabilities forward ================================
  # step_range : (2 * max_dist + 1)^2 
  # n_obs      : Number of GPS observations
  # steps      : Number of simulated steps
  current <- array(0, dim = c(step_range, n_obs, sim_steps))
  
  # center     : Center of step_range (center cell of (2 * buffer + 1)
  # Set to 1 at step #1 for each observation because that's where it actually is
  center <- step_range / 2 + 0.5
  current[center, , 1] <- 1
  
  current2 <- current  # DEBUG
  dest2 <- dest

  for (j in 1:(sim_steps - 1)) {
    # Probabilities across entire nbhd for step j
    step_prob <- as.vector(current[, , j]) * attract[]
    step_prob2 <- as.vector(current2[, , j]) * attract2[] # DEBUG
    # dest has same dimensions as nbhd
    # step_prob is a vector of length step_range
    # to_dest is a matrix with the same dimensions as nbhd
    #   rows are neighborhood cells, columns are neighbors of those cells
    #   each entry is the index of the neighbor cell
    #   e.g. if to_dest[1, 2] = 3, then the 2nd neighbor of the 1st cell of
    #   nbhd is the 3rd cell of nbhd.
    dest[] <- step_prob[as.vector(to_dest)]
    dest2[] <- step_prob2[as.vector(to_dest)] # DEBUG

    # Summing probabilities up to step j to generate step j+1
    current[, , j + 1] <- rowSums(dest, na.rm = TRUE)
    current2[, , j + 1] <- rowSums(dest2, na.rm = TRUE) # DEBUG
  }
  
  #### Run likelihood function from fitting with true parameters

  # Calculate log likelihood 
  predictions <- matrix(0, nrow = sim_steps, ncol = n_obs)
  for (i in 1:n_obs) {
    predictions[, i] <- current[obs[i], i, ]
    # returns the probability for the row associated with the next 
    # observation location, for that observation i, across all time steps
  }

  log_likelihood <- rowSums(log(predictions), na.rm = TRUE)
  # log of product is sum of logs

  saveRDS(predictions, paste0("data/output/simulations/p", current_jag, ".RDS"))
  current <- list(current, current2) # DEBUG
  saveRDS(current, paste0("data/output/simulations/current", current_jag, ".RDS")) # DEBUG
  return(-max(log_likelihood, na.rm = TRUE))
  # Return negative of the maximum log likelihood because we want to minimize
  # Lower negative log likelihood = higher likelihood 

  return(list(current, current2)) # DEBUG
}

# Output analysis --------------------------------------------------------------

## Load parameters and likelihood
load_output <- function(name) {
    dir <- paste0("data/output/", name)
    # ll_files <- list.files(dir)[grep("likelihood_", list.files(dir))]
    # par_files <- list.files(dir)[grep("par_out_", list.files(dir))]
    ll_files <- paste0("likelihood_", 1:njag, ".RDS")
    par_files <- paste0("par_out_", 1:njag, ".RDS")
    ll <- load_if_exists(ll_files, dir)
    par <- load_if_exists(par_files, dir)
    return(list(unlist(ll), par))
}

par_to_df <- function(par) {
    df <- do.call(rbind, lapply(par, function(x) {
        print(x[[1]])
    }))
}

# Simulation -------------------------------------------------------------------

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
    outdf <- rast_df(out)
    return(list(raster = out, df = outdf))
}

# Generate a jaguar path of n steps starting from (x0, y0) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x0, y0, nstep, par, neighb, type = 2, tprob) {
    # type: 1 = env1 only, 2 = multi-state model
    # if (!(x0 %in% 1:100) || !(y0 %in% 1:100)) {
    #     print("Jaguar out of bounds")
    #     return(NULL)
    # }
    path <- matrix(NA, nrow = nstep, ncol = 6)
    state0 <- 1 # Beginning state, irrelevant if type = 1

    path[1, ] <- c(x0, y0, NA, state0, NA, NA)
    
    if (type == 2) {
        # Transition probabilities: p12, p21, p11, p22
        tprob <- c(tprob, 1 - tprob)
    }

    for (i in 2:nstep) {
        pos <- path[i - 1, 1:2]
        state <- path[i - 1, 4]
        if (type == 2 && i %% sim_interval == 0) {
            # Transition to new state
            if (state == 1) {
                if (runif(1) < tprob[1]) state <- 2
            } else {
                if (runif(1) < tprob[2]) state <- 1
            }
        }
        nbhd <- as.vector(make_nbhd(r = env01[[1]], rdf = env01[[2]], sz = neighb,
                          i = cellFromRowCol(env01[[1]], pos[1], pos[2])))
        a1 <- exp(env01[[1]][nbhd] * par[1])$sim1
        if (any(is.na(a1))) a1[is.na(a1)] <- 0
        a1 <- a1 / sum(a1)
        a2 <- exp(env02[[1]][nbhd] * par[2])$sim1
        if (any(is.na(a2))) a2[is.na(a2)] <- 0
        a2 <- a2 / sum(a2)
        
        # par(mfrow = c(1, 2))
        # plot(a1)
        # plot(a2)
        # readLines()

        attract <- switch(state, a1, a2)
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- c(rowColFromCell(env01[[1]], nbhd[step]), 
                       nbhd[step], state, a1[step], a2[step])
    }
    path <- as.data.frame(path)
    names(path) <- c("x", "y", "cell", "state", "a1", "a2")
    return(path)
}

vgram <- function(path, cut = 100) {
    var <- sapply(1:cut, function(t) {
        p1 <- path[1:(nrow(path) - t), 1:2]
        p2 <- path[(t + 1):nrow(path), 1:2]

        out <- sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2)
        return(mean(out))
    })
    return(var)
}

# Plot landscape r with jaguar path and vgram
plot_path <- function(path, vgram = FALSE, new = TRUE, ...) {
    # par(mfrow = c(1, ifelse(vgram, 2, 1)))
    # par(mfrow = c(1, 2))
    path <- path[seq(1, nrow(path), sim_interval), ]

    col1 <- rgb(1, 0, 0, .5)
    col2 <- rgb(0, 0, 1, .8)
    # Plotting environmental variables + path
    if (new) terra::plot(env01[[1]])
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
    # points(path, col = "red", pch = 19, cex = 0.5)
    # lines(path, col = "red")

    # Plotting variogram
    # if (!vgram) return(NULL)
    # plot(vgram(path, ...), type = "l", xlab = "Time lag", ylab = "Variance")
}

# Data =========================================================================

# WGS84 projection
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Jaguar movement data, ID numbers, and metadata
jag_move <- readRDS("data/jag_data_BR.RDS")
jag_id <- readRDS("data/jag_list.RDS"); njag <- nrow(jag_id)
jag_meta <- data.table(read.csv("data/input/jaguars/jaguar_metadata.csv"))

# RasterStack of environmental variables 
# see 01_generate_data.R for details
brazil_ras <- rast("data/env_layers.grd")
# RasterStack of environmental variables, but as a data frame ('brdf')
load("data/env_layers.RData")

msg("Loaded data")
