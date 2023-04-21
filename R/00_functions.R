## Basic functions and libraries

# Libraries
library(raster) # update to terra at some point
library(tidyverse)
library(data.table)
# library(ctmm)
# library(finch)
# library(amt) 
# library(mixtools)
# library(gstat)

# Functions ====================================================================

# General ----------------------------------------------------------------------

# Save raster x as filename fn under ./data/
save_ras <- function(x, fn) {
    # Save raster often because it can get corrupted while working
    writeRaster(x, paste0("data/", fn), format = "raster",
                overwrite = TRUE)
}

# Output message m in logfile f
msg <- function(m, f = "data/output/run_log.txt") {
    m <- paste(format(Sys.time(), "%d.%m.%y  %R"), m)
    print(m)
    cat(m, file = f, append = TRUE, sep = "\n")
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

# Movement model ---------------------------------------------------------------

# Normalize probabilities across neighbors of each cell
norm_nbhd <- function(v) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = T)
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
  attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
                   par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  attract_h <- exp(par[7] * env$home)
  attract <- norm_nbhd(attract_e) * norm_nbhd(attract_h) #* norm_nbhd(attract_t)

  # Attractiveness function 2: just home range ---------------------------------
  # attract_h <- exp(par[1] * env$home)
  # attract <- norm_nbhd(attract_h) 

  # Attractiveness function 3: turn angle --------------------------------------
  # attract_t <- exp(par[8] * turn) # think about functional form of h & t


  # Array for propagating probabilities forward ================================
  # step_range : (2 * buffersize + 1)^2 (= 529)
  # n_obs      : Number of GPS observations
  # steps      : Number of simulated steps
  current <- array(0, dim = c(step_range, n_obs, steps))

  # center     : Center of step_range (center cell of (2 * buffer + 1)
  # Set to 1 at step #1 for each observation because that's where it actually is
  center <- step_range / 2 + 0.5
  current[center, , 1] <- 1
  for (j in 1:(steps - 1)) {
    # Probabilities across entire nbhd for step j
    step_prob <- as.vector(current[, , j]) * attract[]

    # dest has same dimensions as nbhd
    # step_prob is a vector of length step_range
    # to_dest is a matrix with the same dimensions as nbhd
    #   rows are neighborhood cells, columns are neighbors of those cells
    #   each entry is the index of the neighbor cell
    #   e.g. if to_dest[1, 2] = 3, then the 2nd neighbor of the 1st cell of
    #   nbhd is the 3rd cell of nbhd.
    dest[] <- step_prob[as.vector(to_dest)]
    # Summing probabilities up to step j to generate step j+1
    current[, , j + 1] <- rowSums(dest, na.rm = T)

  }
  
  # Calculate log likelihood 
  predictions <- matrix(0, nrow = steps, ncol = n_obs)
  for (i in 1:n_obs) {
    predictions[, i] <- current[obs[i], i, ]
    # returns the probability for the row associated with the next 
    # observation location, for that observation i, across all time steps
  }
  log_likelihood <- rowSums(log(predictions), na.rm = T)
  # log of product is sum of logs
  return(-max(log_likelihood))
  # Return negative of the maximum log likelihood because we want to minimize
  # Lower negative log likelihood = higher likelihood 
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
      r, rdf$row[i] + ind1[j],
      rdf$col[i] + ind2[j]
    )
  }
  return(mat)
}

# Generates a random field with a given correlation structure
# b: beta parameter for gstat, s: sill, r: range, n: nugget
gen_landscape <- function(size = 100, b = 1, s = 0.03, r = 10, n = 0) {   
    xy <- expand.grid(1:size, 1:size); names(xy) <- c("x", "y")

    # Autocorrelation model
    model <- gstat(formula = z ~ 1, locations = ~x + y, dummy = T, beta = b, 
                   model = vgm(psill = s, range = r, model = "Exp", nugget = n), 
                   nmax = 20)
    out <- predict(model, newdata = xy, nsim = 1)
    if (any(out < 0)) out[out < 0] <- 0

    # Output as both raster and data frame
    gridded(out) <- ~x + y; out <- raster(out)
    raster::plot(out)
    outdf <- as.data.frame(out)
    outdf <- cbind(outdf, rowColFromCell(out, seq_len(nrow(outdf))))
    return(list(raster = out, df = outdf))
}

# Generate a jaguar path of n steps starting from (x0, y0) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x0, y0, nstep, par = c(1, 1), neighb = 5) {
    if (!(x0 %in% 1:100) | !(y0 %in% 1:100)) {
        print("Jaguar out of bounds")
        return(NULL)
    }
    path <- matrix(NA, nrow = nstep, ncol = 2)
    x <- x0; y <- y0; path[1, ] <- c(x, y)
    for (i in 2:nstep) {
        pos <- path[i - 1, ]
        nbhd <- make_nbhd(r = env1[[1]], rdf = env1[[2]], sz = neighb,
                          i = cellFromRowCol(env1[[1]], pos[1], pos[2]))
        attract <- exp(env1[[1]][nbhd] * par[1]) # + env2[[1]][nbhd] * par[2]
        if (any(is.na(attract))) attract[is.na(attract)] <- 0
        attract <- attract / sum(attract)
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- rowColFromCell(env1[[1]], nbhd[step])
    }
    path <- as.data.frame(path)
    names(path) <- c("x", "y")
    return(path)
}

vgram <- function(path, cut = 100) {
    var <- sapply(1:cut, function(t) {
        p1 <- path[1:(nrow(path) - t),]
        p2 <- path[(t + 1):nrow(path),]

        out <- sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2)
        return(mean(out))
    })
    return(var)
}

# Plot landscape r with jaguar path and vgram
plot_path <- function(path, ...) {
    par(mfrow = c(1, 2))

    # Plotting environmental variables + path
    raster::plot(env1[[1]])
    points(path, col = "red", pch = 19, cex = 0.5)
    lines(path, col = "red")
    # raster::plot(env2[[1]])
    # points(path, col = "red", pch = 19, cex = 0.5)
    # lines(path, col = "red")

    # Plotting variogram
    plot(vgram(path, ...), type = "l", xlab = "Time lag", ylab = "Variance")
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
brazil_ras <- stack("data/env_layers.grd")
# RasterStack of environmental variables, but as a data frame ('brdf')
load("data/env_layers.RData")

msg("Loaded data")

# ==============================================================================