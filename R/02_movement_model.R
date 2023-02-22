### Might need to add supercomputer functionality back in.
## !IMPORTANT!: Run this through master.R

## Functions ===================================================================

# Outputs a matrix of cell numbers corresponding to brazil_ras
# based on a central cell (i) and a buffer size around that cell (sz)
make_nbhd <- function(i, sz) {
  # if x is the center square and o are neighbors, and e.g. sz = 2
  # (2*sz + 1)^2 represents the following:
  # o o o o o
  # o o o o o
  # o o x o o
  # o o o o o
  # o o o o o

  mat <- matrix(0, nrow = length(i), ncol = (2 * sz + 1)^2)

  # values to add to central cell's row/col to get neighborhood cells' row/col
  ind1 <- t(rep(-sz:sz, each = 2 * sz + 1))
  ind2 <- t(rep(-sz:sz, 2 * sz + 1))
  for (j in seq_len(length(ind1))) {
    mat[, j] <- cellFromRowCol(
      brazil_ras, brdf$row[i] + ind1[j],
      brdf$col[i] + ind2[j]
    )
  }
  return(mat)
}

# For vector v with entries corresponding
norm_nbhd <- function(v) {
  out <- matrix(v[nbhd], nrow = nrow(nbhd), ncol = ncol(nbhd))
  out <- out / rowSums(out, na.rm = T)
}

# Returns negative of the maximum log likelihood given a set of parameters (par)  
loglike_fun <- function(par) {
  # par        : Initial values of parameters for optimization
  # nbhd       : Neighborhood
  # step_range : Step range
  # n_obs      : Number of 
  # steps:     : Number of steps simulated
  # to_dest    : For each cell of the extended neighborhood of the path, what 
  #               are the immediate neighbors? Rows are path cells, columns are 
  #               neighbors.
  # obs        :  

  # Utility function for environmental variables
  # utility <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #   par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])

  # Attractiveness function 
  attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
                   par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  attract_h <- exp(par[7] * home)
  attract_t <- exp(par[8] * turn)

  # 'attract' here is basically the pull factor for each cell of nbhd
  attract <- norm_nbhd(attract_e) * norm_nbhd(attract_h) * norm_nbhd(attract_t)

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

    # basically want to do tmp1 = tmp[to_dest]
    dest[] <- step_prob[as.vector(to_dest)]
    current[, , j + 1] <- rowSums(dest, na.rm = T)
    # Summing probabilities up to step j to generate step j+1
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

ncell <- (buffersize * 2 + 1)^2
msg(paste("Making", ncell, "cell neighborhood for each cell in Brazil"))
neigh <- make_nbhd(seq_len(nrow(brdf)), buffersize)

### Fitting home ranges 
if (refit_homes) {
  for (i in jag_id) {
    jag_traject <- jag_move[ID == as.numeric(i)]
    jag_traject <- as.telemetry(jag_traject, timeformat = "auto")

    # Calculate and plot variogram for individual
    var_jag_traject <- variogram(jag_traject); plot(var_jag_traject)
    # Guess a model for individual
    guess <- ctmm.guess(jag_traject, interactive = F)
    # Fit data to best-guess model
    jag_traject_guess <- ctmm.fit(jag_traject, guess)
    # Fit and plot autocorrelated kernel density estimate
    jag_kde <- akde(jag_traject, jag_traject_guess)
    plot(jag_kde); plot(jag_traject, add = TRUE)
  }
}

### Fitting turn angle distributions
if (refit_turns) {
  for (i in jag_id) {
    jag_traject <- jag_move[ID == as.numeric(i), 3:4]

    ## Fit mixture model here

  }
}


### Fitting

msg("Fitting each jaguar separately")
n_iter <- length(jag_id) # Number of iterations (i.e. number of jaguars)

for (i in i_initial:n_iter) {

  msg(paste0("Jaguar #: ", i, " / ", n_iter))
  # Observed trajectory of jaguar i
  jag_traject <- jag_move[ID == as.numeric(jag_id[i]), 3:4]
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
  nbhd_index <- make_nbhd(jag_traject_cells, size_out)
  # Each entry in the list is the immediate neighborhood of each cell in the 
  # extended neighborhood, as represented by a cell number of brazil_ras
  nbhd_list <- lapply(seq_len(nrow(nbhd_index)), function(x) {
    cells <- nbhd_index[x, ]
    out <- neigh[cells, ]; row.names(out) <- as.character(cells)
    return(out)
  })
  nbhd <- do.call(rbind, nbhd_list)
  # Reindexing allows linkage of row numbers from nbhd to brazil_ras cells
  nbhd_index <- as.vector(t(nbhd_index))

  msg("Getting indices of immediate neighborhood of each cell...")
  # For each cell of the extended neighborhood of the path, what are
  # the immediate neighbors? Rows are path cells, columns are neighbors.
  # All row lengths standardized by turning missing neighbors into NAs.
  to_dest <- tapply(seq_len(length(nbhd)), nbhd, function(x) {
    out <- c(x, rep(NA, ncol(nbhd) - length(x)))
    return(out)
  })
  to_dest <- t(matrix(unlist(to_dest), nrow = ncol(nbhd), ncol = nrow(nbhd)))
  dest <- matrix(0, nrow = nrow(nbhd), ncol = ncol(nbhd))

  # Normalizing desired environmental variables for extended neighborhood
  env <- brdf[nbhd_index, ]
  env <- sweep(env, 2, colMeans(env), "-")       # Subtract mean
  env <- sweep(env, 2, apply(env, 2, sd), "/")   # Divide by std dev
  row.names(env) <- seq_len(length(nbhd_index))  # same indexing as env

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
    # msg(which(index_mat[,y]==test[num]))
    obs[y] <- which(index_mat[, y] == test[num])
    # what is this actually doing I still don't know
  } 

  step_range <- nrow(nbhd) / length(jag_traject_cells)
  n_obs <- length(jag_traject_cells)
  steps <- 25

  if (model_calcnull) {
    msg(paste0("Calculating null likelihood for jaguar ", i))
    null_likelihood <- loglike_fun(c(rep(0, 6)))
    saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
  }

  if (model_fiteach) {
    param <- rnorm(6)
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

          cat(paste0("jaguar ", i, " fitted ", date()),
            file = "data/output/run_log.txt",
            append = TRUE, sep = "\n"
          )
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
