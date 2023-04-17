### Might need to add supercomputer functionality back in.
## !IMPORTANT!: Run 00_basics.R first to set up the environment

## Switches ====================================================================

refit_homes     <- F              # Refit home ranges (AKDE) 
refit_turns     <- F              # Refit turn distributions (MM)
refit_model     <- T              # Refit movement model parameters
model_calcnull  <- T              # Calculate null likelihoods 
                                    # refit_model must be TRUE for this one

npar            <- 7              # Number of parameters in current model

i_initial       <- 1              # Individual to start at
buffersize      <- 1              # Jaguars move 1px (1km) at a time
n_iter          <- nrow(jag_id)   # Number of individuals

# Functions ====================================================================

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

  # Attractiveness function 1: environmental variables + home range
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  # attract_h <- exp(par[7] * env$home)
  # attract <- norm_nbhd(attract_e) * norm_nbhd(attract_h) # * norm_nbhd(attract_t)

  # Attractiveness function 2: just home range
  attract_h <- exp(par[1] * env$home)
  attract <- norm_nbhd(attract_h) 

  # Attractiveness function 3: turn angle
  # attract_t <- exp(par[8] * turn) # think about functional form of h & t

  # Array for propagating probabilities forward
  # step_range : (2 * buffersize + 1)^2 (= 9)
  # n_obs      : Number of observations
  # steps      : Number of simulated steps
  current <- array(0, dim = c(step_range, n_obs, steps))

  # center     : Center of step_range (center cell of (2 * buffer + 1)
  # Set to 1 at step #1 for each observation because that's where it actually is
  center <- step_range / 2 + 0.5
  current[center, , 1] <- 1
  for (j in 1:(steps - 1)) {
    # Probabilities across all step range for step j
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
  # current has everything - but need to know where the next obs was (row) for
  # each column, and will sum across each time steps
  predictions <- matrix(0, nrow = steps, ncol = n_obs)
  for (i in 1:n_obs) {
    predictions[, i] <- current[obs[i], i, ]
    # returns back the probability for the row
    # associated with the next observation location, for that observation i,
    # across all time steps
  }
  log_likelihood <- rowSums(log(predictions), na.rm = T)
  return(-max(log_likelihood))
}

### Body =======================================================================

### Fitting home ranges
if (refit_homes) {
  ## Potentially: run the main fitting algorithm with and without home ranges
  ## and use model selection to decide whether to include home ranges on a per
  ## individual basis (?) 
  for (id in jag_id$jag_id) {
    msg(paste0("Fitting home range for individual ", id))
    # Individual trajectory
    jag_traject <- as.telemetry(jag_move[ID == as.numeric(id)],
                                timeformat = "auto")

    # Calculate and plot variogram for individual
    var_jag_traject <- ctmm::variogram(jag_traject)
    # Guess a model for individual
    guess <- ctmm.guess(jag_traject, interactive = F)
    # Fit data to best-guess model
    jag_traject_guess <- ctmm.fit(jag_traject, guess)
    # Fit and plot autocorrelated kernel density estimate
    jag_kde <- raster(akde(jag_traject, jag_traject_guess))
    jag_kde <- projectRaster(jag_kde, brazil_ras)
    # Saving plots
    # name <- paste0("data/output/homeranges/homerange_", id, ".png")
    # png(filename = name, width = 1000, height = 500)
    # par(mfrow = c(1, 2))
    # plot(var_jag_traject)
    # plot(jag_kde); plot(jag_traject, add = TRUE)   
    # dev.off()
    jag_kde <- resample(jag_kde, brazil_ras[[1]], method = "ngb")
    save_ras(jag_kde, paste0("output/homeranges/homerange_", id, ".tif"))

  }
}

### Fitting turn angle distributions 
if (refit_turns) {
  msg("Fitting mixture model for turn angle distributions")
  turn_models <- lapply(jag_id$jag_id, function(id) {
    msg(paste("Fitting turn angle mixture model for individual", id))
    # Individual trajectory (using amt::make_track for turn angles)
    jag_traject <- as.telemetry(jag_move[ID == as.numeric(id)],
                                timeformat = "auto")
    jag_track <- make_track(jag_traject, x, y)
    turns <- direction_abs(jag_track)

    # Fitting a normal mixture model
    mm_turns <- normalmixEM(turns[-which(is.na(turns))])

    name <- paste0("data/output/turnmodels/turnmodel_", id, ".png")
    png(filename = name, width = 500, height = 500)
    plot.mixEM(mm_turns, whichplots = 2)
    dev.off()
  })
  # saveRDS(turn_models, "data/output/turn_models.RDS")  
}

### Fitting environmental parameters
if (refit_model) {
  msg("Fitting model parameters")
  ncell <- (buffersize * 2 + 1)^2
  msg(paste("Making", ncell, "cell neighborhood for each cell in Brazil"))
  nbhd0 <- make_nbhd(i = seq_len(nrow(brdf)), sz = buffersize)                   # 6.4s

  for (i in i_initial:n_iter) {

    msg(paste0("Jaguar #: ", i, " / ", n_iter))
    # Observed trajectory of jaguar i
    id <- as.numeric(jag_id[i])
    jag_traject <- jag_move[ID == id, 3:4]
    jag_traject_cells <- cellFromXY(brazil_ras, jag_traject)
    # Calculating step distances; divide by cell size then take hypotenuse
    dist <- (jag_traject[-nrow(jag_traject), ] - jag_traject[-1, ]) /
            xres(brazil_ras)
    dist <- (rowSums(dist^2))^.5
    # Making neighborhoods for each point in trajectory, with buffer size as 
    # twice maximum step distance
    max_dist <- max(dist)             
    size_out <- ceiling(max_dist * 2) 

    msg("Building neighborhoods for each cell")
    # Extended neighborhoods of each cell in individual's trajectory
    nbhd_index <- make_nbhd(i = jag_traject_cells, sz = size_out)
    # Each entry in the list is the immediate neighborhood of each cell in the 
    # extended neighborhood, as represented by a cell number of brazil_ras
    nbhd_list <- lapply(seq_len(nrow(nbhd_index)), function(i) {                 # 13.5s
      row_num <- seq_len(ncol(nbhd_index)) + (i - 1) * ncol(nbhd_index)
      names(row_num) <- nbhd_index[i, ] # index names for i
      out <- matrix(row_num[as.character(nbhd0[nbhd_index[i, ], ])], 
                    nrow = length(row_num), ncol = ncol(nbhd0))
      return(out)
    })
    nbhd <- do.call(rbind, nbhd_list)
    # Reindexing allows linkage of row numbers from nbhd to brazil_ras cells
    nbhd_index <- as.vector(t(nbhd_index))

    msg("Getting indices of immediate neighborhood of each cell...")
    # For each cell of the extended neighborhood of the path, what are
    # the immediate neighbors? Rows are path cells, columns are neighbors.
    # All row lengths standardized by turning missing neighbors into NAs.
    to_dest <- tapply(seq_len(length(nbhd)), nbhd, function(x) {                 # 36s
      # print(x)
      out <- c(x, rep(NA, ncol(nbhd) - length(x)))
      return(out)
    })
    to_dest <- t(matrix(unlist(to_dest), nrow = ncol(nbhd), ncol = nrow(nbhd)))
    dest <- matrix(0, nrow = nrow(nbhd), ncol = ncol(nbhd))

    # Adding individual home range (AKDE) to brdf
    home <- raster(paste0("data/output/homeranges/homerange_", id, ".grd"))
    brdf$home <- as.vector(home)

    # Normalizing desired environmental variables for extended neighborhood
    env <- brdf[nbhd_index, ]
    env[, 1:6] <- sweep(env[, 1:6], 2, colMeans(env)[1:6], "-") 
    env[, 1:6] <- sweep(env[, 1:6], 2, apply(env, 2, sd)[1:6], "/") 
    # Make indexing consistent with env
    row.names(env) <- seq_len(length(nbhd_index))

    # Building observed data to test against
    index_mat <- matrix(
      data = seq_len(length(nbhd_index)),
      nrow = (nrow(nbhd) / length(jag_traject_cells)),
      ncol = length(jag_traject_cells)
    )
    obs <- vector(length = ncol(index_mat) - 1)
    for (y in 1:(ncol(index_mat) - 1)) {
      test <- which(nbhd_index == jag_traject_cells[y + 1])
      num <- which(index_mat[1, y] < test & test < index_mat[nrow(index_mat), y])
      obs[y] <- which(index_mat[, y] == test[num])
    } 

    step_range <- nrow(nbhd) / length(jag_traject_cells)
    n_obs <- length(jag_traject_cells)
    steps <- 25

    # Calculate null likelihoods for each jaguar if not already done
    if (model_calcnull) {
      msg(paste0("Calculating null likelihood for jaguar ", i))
      null_likelihood <- loglike_fun(c(rep(0, npar)))
      saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
    } else {
      param <- rnorm(npar)
      msg("Running optim...")
      ntries <- 0
      ## Main fitting loop, tries each individual 20 times and moves on if no fit
      while (ntries <= 20) {
        tryCatch({
            par_out <- optim(param, loglike_fun)
            saveRDS(par_out, paste0("data/output/par_out_", i, ".RDS"))

            msg("Running loglike_fun...")
            likelihood <- loglike_fun(par_out[[1]])
            saveRDS(likelihood, paste0("data/output/likelihood_", i, ".RDS"))

            msg(paste0("jaguar ", i, " fitted ", date()))
            ntries <- 21 # End while loop
          },
          error = function(e) {
            msg(e)
            msg(paste("Try #:", ntries))
            if (ntries == 20) {
              msg("Skipping, couldn't fit in 20 tries")
            } else {
              msg("Retrying")
            }
          },
          finally = {
            ntries <- ntries + 1
          }
        )
      }
    }
  }
}