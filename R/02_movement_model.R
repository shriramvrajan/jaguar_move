### Might need to add supercomputer functionality back in.
## !IMPORTANT!: Run 00_basics.R first to set up the environment

## Switches ====================================================================

refit_homes     <- FALSE            # Refit home ranges (AKDE) 
refit_turns     <- FALSE            # Refit turn distributions (MM)
refit_model     <- TRUE             # Refit movement model parameters
model_calcnull  <- FALSE            # Calculate null likelihoods 
                                    # refit_model must be TRUE for this one

npar            <- 7              # Number of parameters in current model
steps           <- 25             # How many steps to simulate forward

i_initial       <- 1              # Individual to start at
buffersize      <- 1              # Jaguars move 1px (1km) at a time
n_iter          <- nrow(jag_id)   # Number of individuals

### Body =======================================================================

### Fitting home ranges
if (refit_homes) {
  ## Potentially: run the main fitting algorithm with and without home ranges
  ## and use model selection to decide whether to include home ranges on a per
  ## individual basis (?) 
  for (id in jag_id$jag_id) {
    message(paste0("Fitting home range for individual ", id))
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
    save_raster(jag_kde, paste0("output/homeranges/homerange_", id, ".tif"))

  }
}

### Fitting turn angle distributions 
if (refit_turns) {
  message("Fitting mixture model for turn angle distributions")
  turn_models <- lapply(jag_id$jag_id, function(id) {
    message(paste("Fitting turn angle mixture model for individual", id))
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
  message("Fitting model parameters")
  ncell <- (buffersize * 2 + 1)^2
  message(paste("Making", ncell, "cell neighborhood for each cell in Brazil"))
  nbhd0 <- make_nbhd(i = seq_len(nrow(brdf)), sz = buffersize)                   # 6.4s

  for (i in i_initial:n_iter) {

    message(paste0("Jaguar #: ", i, " / ", n_iter))
    id <- as.numeric(jag_id[i])

    # Adding individual home range (AKDE) to brdf
    home <- rast(paste0("data/homeranges/homerange_", id, ".grd"))
    brdf$home <- as.vector(home)
    # Observed trajectory of jaguar i
    jag_traject <- jag_move[ID == id, 3:4]
    jag_traject_cells <- cellFromXY(brazil_ras, jag_traject)
    n_obs <- length(jag_traject_cells)
    # Calculating step distances; divide by cell size then take hypotenuse
    dist <- (jag_traject[-nrow(jag_traject), ] - jag_traject[-1, ]) /
            xres(brazil_ras)
    dist <- (rowSums(dist^2))^.5
    # Making neighborhoods for each point in trajectory, with buffer size as 
    # twice maximum step distance
    max_dist <- ceiling(max(dist) * 2)           
    step_range <- (max_dist * 2 + 1) ^ 2

    # Prepare input for fitting, given trajectory, neighborhood size, and number
    # of steps to simulate forward
    prep_model_objects(jag_traject_cells, max_dist, nsteps, r = brazil_ras, rdf = brdf,
               nbhd0 = nbhd0)

    # Normalizing desired environmental variables for extended neighborhood
    envdf <- brdf[, c(1:6, 10)]
    env <- envdf[nbhd_index, ]
    env <- sweep(env, 2, colMeans(env), "-") 
    env <- sweep(env, 2, apply(env, 2, sd), "/") 
    # Make indexing consistent with env
    row.names(env) <- seq_len(length(nbhd_index))

    # Calculate null likelihoods for each jaguar if not already done
    if (model_calcnull) {
      message(paste0("Calculating null likelihood for jaguar ", i))
      null_likelihood <- log_likelihood(c(rep(0, npar)))
      saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
    } else {
      param <- rnorm(npar)
      message("Running optim...")
      ntries <- 0
      ## Main fitting loop, tries each individual 20x and moves on if no fit
      while (ntries <= 20) {
        tryCatch({
            par_out <- optim(param, loglike_fun)
            saveRDS(par_out, paste0("data/output/par_out_", i, ".RDS"))

            message("Running loglike_fun...")
            likelihood <- log_likelihood(par_out[[1]])
            saveRDS(likelihood, paste0("data/output/likelihood_", i, ".RDS"))

            message(paste0("jaguar ", i, " fitted ", date()))
            ntries <- 21 # End while loop
          },
          error = function(e) {
            message(e)
            message(paste("Try #:", ntries))
            if (ntries == 20) {
              message("Skipping, couldn't fit in 20 tries")
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
  }
}