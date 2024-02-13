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
    jag_track <- make_track(id)
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
  ncell <- (buffersize * 2 + 1)^2
  message(paste("Making", ncell, "cell neighborhood for each cell in Brazil"))
  nbhd0 <- make_nbhd(i = seq_len(nrow(brdf)), sz = buffersize)                   # 6.4s

  # foreach(i = i_todo) %dopar% {
  for (i in i_todo) {
    message(paste0("Jaguar #: ", i))
    id <- as.numeric(jag_id[i])

    # Adding individual home range (AKDE) to brdf
    home <- rast(paste0("data/homeranges/homerange_", id, ".grd"))
    brdf$home <- as.vector(home)
    envdf <- brdf[, c(1:6, 10)]

    # Observed trajectory of jaguar i
    jag_traject <- jag_move[ID == id, 3:4]
    
    if (holdout_set && nrow(jag_traject) > 100) {
      hold <- seq_len(ceiling(nrow(jag_traject) * holdout_frac))
      jag_traject <- jag_traject[hold, ]
    }

    jag_traject_cells <- cellFromXY(brazil_ras, jag_traject)
    n_obs <- length(jag_traject_cells)
    # Calculating step distances; divide by cell size then take hypotenuse
    dist <- (jag_traject[-nrow(jag_traject), ] - jag_traject[-1, ]) /
            xres(brazil_ras)
    dist <- (rowSums(dist^2))^.5
    
    if (refit_model0 == TRUE) {
      # Traditional SSF for comparison
      max_dist <- floor(max(dist)) 
      nbhd <- make_nbhd(i = jag_traject_cells, sz = max_dist)
      obs <- unlist(lapply(seq_len(nrow(nbhd) - 1), function(step) {
        which(nbhd[step, ] == jag_traject_cells[step + 1])
      }))
      nbhd_index <- as.vector(t(nbhd))
      track <- make_track(id)
      sl_emp <- as.vector(na.exclude(track$sl))
      ta_emp <- as.vector(na.exclude(track$ta))
      mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, max_dist = max_dist)
    } else {
      # New model
      max_dist <- ceiling(max(dist) * 2)   
      prep_model_objects(jag_traject_cells, max_dist, r = brazil_ras, rdf = brdf, 
                        nbhd0 = nbhd0)
    }
  
    # Normalizing desired environmental variables for extended neighborhood
    env <- envdf[nbhd_index, ]
    env <- sweep(env, 2, colMeans(env), "-") 
    env <- sweep(env, 2, apply(env, 2, sd), "/") 
    # Make indexing consistent with env
    row.names(env) <- seq_len(length(nbhd_index))

    # Model objects as list
    objects <- ifelse(refit_model0, 
                      list(env, max_dist, mk, obs),
                      list(env, nbhd, max_dist, sim_steps, to_dest, obs))
    # browser()

    # Calculate null likelihoods for each jaguar if not already done
    if (model_calcnull) {
      message(paste0("Calculating null likelihood for jaguar ", i))
      null_likelihood <- loglike(c(rep(0, npar)), objects)
      saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
    } else {
      param <- rnorm(npar)
      message("Running optim...")
      run_optim(param, objects, i)
    } 
  }
}
