### Might need to add supercomputer functionality back in.
## !IMPORTANT!: Run through master.R 

### Fitting home ranges --------------------------------------------------------
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
    guess <- ctmm.guess(jag_traject, interactive = FALSE)
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

### Fitting turn angle distributions -------------------------------------------
if (refit_turns) {
  message("Fitting mixture model for turn angle distributions")
  turn_models <- lapply(jag_id$jag_id, function(id) {
    message(paste("Fitting turn angle mixture model for individual", id))
    # Individual trajectory 
    jag_track <- make_full_track(id)
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

### Fitting environmental parameters -------------------------------------------
if (refit_model) {

  message("Fitting model parameters")

  foreach(i = i_todo) %dopar% {
  # for (i in i_todo) {  # easier to debug
    message(paste0("Jaguar #: ", i))
    id                <- as.numeric(jag_id[i])
    jag_traject       <- jag_move[ID == id, 3:4]
    jag_traject_cells <- cellFromXY(brazil_ras, jag_traject)
    n_obs             <- length(jag_traject_cells)
    # Calculating step distances; divide by cell size then take hypotenuse
    dist     <- (jag_traject[-nrow(jag_traject), ] - jag_traject[-1, ]) /
                 xres(brazil_ras)
    dist     <- (rowSums(dist^2))^.5
    max_dist <- ceiling(max(dist) * 1.5)
    # home      <- rast(paste0("data/homeranges/homerange_", id, ".grd"))
    # brdf$home <- as.vector(home)
    envdf    <- brdf[, c(1:6)] # add 10 for homerange
    
    if (holdout_set && nrow(jag_traject) > 100) {
      hold <- seq_len(ceiling(nrow(jag_traject) * holdout_frac))
      jag_traject <- jag_traject[hold, ]
    }

    param0 <- rnorm(npar)

    if (model_type == 1) {
      message("Using traditional step selection function model")
      nbhd <- make_nbhd(i = jag_traject_cells, sz = max_dist)
      obs <- sapply(seq_along(jag_traject_cells), function(i) {
        if (i == length(jag_traject_cells)) {
          return(NULL)
        } else {
          step <- jag_traject_cells[i + 1]
          return(which(nbhd[i, ] == step))
        }
      }) %>% unlist()
      env <- scale(envdf[unique(nbhd), ])
      if (any(is.na(env))) env[which(is.na(env))] <- 0
      nbhd_c <- matrix(as.character(nbhd), nrow = nrow(nbhd), ncol = ncol(nbhd)) # needs to be character for this one
      objects <- list(nbhd_c, obs, env, max_dist)
    } else if (model_type == 2) {
      message("Using path propagation model")
      objects <- prep_model_objects(jag_traject_cells, max_dist, envdf)
    } else {
      stop("Invalid model type")
    }
    # Calculate null likelihoods for each jaguar if not already done
    if (model_calcnull) {
      message(paste0("Calculating null likelihood for jaguar ", i))
      null_likelihood <- switch(model_type,
                                log_likelihood0(c(rep(0, npar)), objects),
                                log_likelihood1(c(rep(0, npar)), objects))
      saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
    } else {
      message("Running optim...")
      run_optim(param0, objects, i)
    } 

  }
}
