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
  
  ll_func <- switch(model_type, log_likelihood0, log_likelihood)

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
    envdf    <- brdf[, c(1:6)] # add 1 for footprint, 10 for homerange
    
    if (holdout_set && nrow(jag_traject) > 100) {
      hold <- seq_len(ceiling(nrow(jag_traject) * holdout_frac))
      jag_traject <- jag_traject[hold, ]
      jag_traject_cells <- jag_traject_cells[hold]
    }

    # Initial parameter values
    
    param0 <- rep(1, npar)

    # param0 <- load_output("pp1", i)$par
    # p_s <- 1 - exp01(param0[7])
    # param0[length(param0)] <- -log((1/p_s - 1) / 8) # test, only works for step_size=1
    
    # param0 <- load_output("pp1_2")[[i]]$par

    # Preparing model objects based on model type; 1 = SSF, 2 = path propagation
    if (model_type == 1) {
      message("Using traditional step selection function model")

      # empirical movement kernel ----------------------------------------------
      track <- make_full_track(id)
      sl_emp <- as.vector(na.exclude(track$sl))
      ta_emp <- as.vector(na.exclude(track$ta))
      mk <- make_movement_kernel(sl_emp, ta_emp, n = 5000, max_dist = max_dist)
      #-------------------------------------------------------------------------   
      
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
      objects <- list(nbhd_c, obs, env, max_dist, mk)
      names(objects) <- c("nbhd", "obs", "env", "max_dist", "mk")
    } else if (model_type == 2) {
      message("Using path propagation model")
      objects <- prep_model_objects(jag_traject_cells, max_dist, envdf)
      message(object.size(objects) %>% format(units = "Mb"))
    } else {
      stop("Invalid model type")
    }

    # Calculate null likelihoods for each jaguar if not already done
    if (model_nofit) {
      message(paste0("Calculating null likelihood for jaguar ", i))
      likelihood <- ifelse(model_calcnull, 
                           ll_func(c(rep(0, npar)), objects),
                           ll_func(param0, objects))
      name <- ifelse(model_calcnull, "ll_null_", "ll_")
      saveRDS(likelihood, paste0("data/output/", name, i, ".rds"))
    } else {
      message("Running optim...")
      run_optim(param0, objects, i)
      if (debug_02) {
        message("Debugging 02_movement_model.R")
        par1 <- readRDS(paste0("data/output/par_out_", i, ".rds"))$par
        fname <- paste0("data/output/debug_", i, ".rds")
        ll_func(par1, objects, debug = TRUE) %>% saveRDS(fname)
      }
    } 
  }
}
