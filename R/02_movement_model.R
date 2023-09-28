### Might need to add supercomputer functionality back in.
## !IMPORTANT!: Run 00_basics.R first to set up the environment

## Switches ====================================================================

refit_homes     <- FALSE            # Refit home ranges (AKDE) 
refit_turns     <- FALSE            # Refit turn distributions (MM)
refit_model     <- TRUE             # Refit movement model parameters
model_calcnull  <- FALSE            # Calculate null likelihoods 
                                    # refit_model must be TRUE for this one
refit_model0    <- TRUE             # Refit traditional SSF model
                                    
npar            <- 7              # Number of parameters in current model
steps           <- 25             # How many steps to simulate forward

i_initial       <- 1              # Individual to start at
buffersize      <- 1              # Jaguars move 1px (1km) at a time
n_iter          <- nrow(jag_id)   # Number of individuals

### Body =======================================================================

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

### Fitting traditional SSF ----------------------------------------------------
if (refit_model0) {

  message("Fitting traditional SSF model")
  for (id in jag_id$jag_id) {
    
    nbhd <- make_nbhd()

  }


  # message("Fitting movement kernels for all individuals")
  # for (id in jag_id$jag_id) {
  #   message(paste0("Fitting movement kernel for individual ", id))
  #   # Individual trajectory
  #   jag_traject <- make_track(id)
    
  #   sl_emp <- as.integer(na.exclude(jag_traject$sl) / 1000)
  #   ta_emp <- na.exclude(jag_traject$ta)

  #   build_movement_kernel <- function() {
  #     # Generate 1000 random step lengths and turn angles
  #     # Calculate where they fall on grid 
  #     # Save as vector
  #     sl_samp <- sample(sl_emp, 10000, replace = TRUE)
  #     ta_samp <- sample(ta_emp, 10000, replace = TRUE)
  #     avail <- data.frame(sl = sl_samp, ta = ta_samp, x = sl_samp * cos(ta_samp),
  #                         y = sl_samp * sin(ta_samp))
  #     avail$xi <- sapply(avail$x, function(x) ifelse(x > 0, floor(x), ceiling(x)))
  #     avail$yi <- sapply(avail$y, function(y) ifelse(y > 0, floor(y), ceiling(y)))


  #     if (any(is.na(avail$x))) avail <- avail[-which(is.na(avail$x)), ]
      
  #   }
    
  #   counts <- tapply(avail$x, list(avail$xi, avail$yi), length, simplify = F)
  #   neighb <- matrix(nrow = nrow(counts) + 2, ncol = ncol(counts) + 2)
  #   neighb[2:(nrow(counts) + 1), 2:(ncol(counts) + 1)] <- counts
  #   neighb[which(is.na(neighb))] <- min(neighb, na.rm = TRUE)
  #   neighb <- neighb / sum(neighb)

  #   par(mfrow = c(1, 2))
  #   plot(avail$x, avail$y)
  #   # plot(avail$xi, avail$yi)
  #   terra::plot(rast(neighb))

  #   # Generate 1000 random step lengths and turn angles
  #   # Calculate where they fall on grid 
  #   # Save as vector
  # }
}

### Fitting environmental parameters -------------------------------------------
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
    envdf <- brdf[, c(1:6, 10)]

    # Observed trajectory of jaguar i
    jag_traject <- jag_move[ID == id, 3:4]
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
      obs <- sapply(seq_len(nrow(nbhd) - 1), function(step) {
        which(nbhd[step, ] == jag_traject_cells[step + 1])
      })
      nbhd_index <- as.vector(t(nbhd))
      track <- make_track(id)
      sl_emp <- as.vector(na.exclude(track$sl))
      ta_emp <- as.vector(na.exclude(track$ta))
      mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000)
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

    # Calculate null likelihoods for each jaguar if not already done
    if (model_calcnull) {
      message(paste0("Calculating null likelihood for jaguar ", i))
      null_likelihood <- loglike(c(rep(0, npar)))
      saveRDS(null_likelihood, paste0("data/output/null_", i, ".RDS"))
    } else {
      param <- rnorm(npar)
      message("Running optim...")
      run_optim(param)
    } 
  }
}