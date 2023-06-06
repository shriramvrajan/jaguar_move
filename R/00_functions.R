## Basic functions and libraries

# Libraries
library(raster) # update to terra at some point
library(tidyverse)
library(data.table)
library(ctmm)
library(amt) 
library(gstat)
library(terra)
# library(mixtools)

# Functions ====================================================================

# Basic ------------------------------------------------------------------------

# Save raster x as filename fn under ./data/
save_ras <- function(x, fn) {
    # Save raster often because it can get corrupted while working
    writeRaster(x, paste0("data/", fn), overwrite = TRUE)
}

# Output message m both in console and in logfile f
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

# Normalize probabilities across neighbors of each cell
norm_nbhd <- function(v) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = TRUE)
}

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
      # msg(paste(x, "to_dest"))
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
      msg(y)
      # browser()
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
  attract <- norm_nbhd(exp(par[1] * env1)) # + exp(par[2] * env2)

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
    current[, , j + 1] <- rowSums(dest, na.rm = TRUE)
  }
  
  #### Run likelihood function from fitting with true parameters

  # Calculate log likelihood 
  predictions <- matrix(0, nrow = steps, ncol = n_obs)
  for (i in 1:n_obs) {
    predictions[, i] <- current[obs[i], i, ]
    # returns the probability for the row associated with the next 
    # observation location, for that observation i, across all time steps
  }
  log_likelihood <- rowSums(log(predictions), na.rm = TRUE)
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
    gridded(out) <- ~x + y
    out <- raster(out)
    # raster::plot(out)
    outdf <- as.data.frame(out)
    outdf <- cbind(outdf, rowColFromCell(out, seq_len(nrow(outdf))))
    return(list(raster = out, df = outdf))
}

# Generate a jaguar path of n steps starting from (x0, y0) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x0, y0, nstep, par, neighb, type = 2, tprob) {
    # type: 1 = env1 only, 2 = multi-state model
    if (!(x0 %in% 1:100) || !(y0 %in% 1:100)) {
        print("Jaguar out of bounds")
        return(NULL)
    }
    path <- matrix(NA, nrow = nstep, ncol = 6)
    state0 <- 1 # Beginning state, irrelevant if type = 1
    x <- x0
    y <- y0
    path[1, ] <- c(x, y, NA, state0, NA, NA)
    
    if (type == 2) {
        # Transition probabilities: p12, p21, p11, p22
        tprob <- c(tprob, 1 - tprob)
    }

    for (i in 2:nstep) {
        pos <- path[i - 1, 1:2]
        state <- path[i - 1, 4]
        if (type == 2) {
            # Transition to new state
            if (state == 1) {
                if (runif(1) < tprob[1]) state <- 2
            } else {
                if (runif(1) < tprob[2]) state <- 1
            }
        }
        nbhd <- make_nbhd(r = env1[[1]], rdf = env1[[2]], sz = neighb,
                          i = cellFromRowCol(env1[[1]], pos[1], pos[2]))
        a1 <- exp(env1[[1]][nbhd] * par[1])
        if (any(is.na(a1))) a1[is.na(a1)] <- 0
        a1 <- a1 / sum(a1)
        a2 <- exp(env2[[1]][nbhd] * par[2])
        if (any(is.na(a2))) a2[is.na(a2)] <- 0
        a2 <- a2 / sum(a2)
        attract <- switch(state, a1, a2)
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- c(rowColFromCell(env1[[1]], nbhd[step]), 
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
plot_path <- function(path, vgram = T, new = T, ...) {
    par(mfrow = c(1, ifelse(vgram, 2, 1)))

    col1 <- rgb(1, 0, 0, .5)
    col2 <- rgb(0, 0, 1, .8)
    # Plotting environmental variables + path
    if (new) raster::plot(env1[[1]])
    points(path, col = c(col1, col2)[path$state], pch = 19, cex = 0.5)
    for (i in 1:(length(path$state) - 1)) {
        segments(path$x[i], path$y[i], path$x[i + 1], path$y[i + 1], 
                 col = c(col1, col2)[path$state[i]])
    }
    # raster::plot(env2[[1]])
    # points(path, col = "red", pch = 19, cex = 0.5)
    # lines(path, col = "red")

    # Plotting variogram
    if (!vgram) return(NULL)
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