rm(list = ls())
source("R/functions.R")     # Existing functions
source("R/classes.R")       # New classes

set.seed(7)

fit_individuals <- TRUE
test_holdout    <- FALSE

# 1 for step selection, 2 for path propagation
model_type <- 2

if (fit_individuals) {
    config <- list(   
        # Model parameters
        model_type        = model_type,
        npar              = switch(model_type, 9, 9),    # Number of parameters
        step_size         = 1,    # Minimum step size (inner neighborhood) in pixels
        n_jump_range      = 0:4,  # Range of jump sizes to consider
        individuals       = 12, # NULL = all individuals, or vector of IDs

        # Holdout set parameters
        holdout_set  = FALSE,     # Whether to reserve holdout set (T/F)
        holdout_frac = 0.7,       # Proportion of data to use for training

        # Parallel processing parameters
        parallel = FALSE,          # Whether to use parallel processing (T/F)
        n_cores  = 6,             # Number of cores to use if parallel

        # Model fitting options
        fit_model      = TRUE,    # Whether to fit the model (T/F)
        model_calcnull = FALSE    # Whether to calculate null model likelihood (T/F)
    )

    batch <- empirical_batch$new(config)
    results <- batch$run_all()

    message("Saving results...")
    saveRDS(results, paste0("data/output/empirical_results_", 
                            switch(config$model_type, "ss", "pp"), "_",
                            Sys.Date(), ".rds"))
}

## Holdout set evaluation ======================================================

if (test_holdout) {
    ss_model <- step_selection_model$new()
    pp_model <- path_propagation_model$new()

    ### Change seq_len(nrow(res)) part
    ll_holdout <- sapply(seq_len(nrow(res)), function(i) {
        id <- res$ID[i]
        message(paste0("Holdout analysis for individual ", id))

        par_ss_i <- as.numeric(res[i, 3:11])
        par_pp_i <- as.numeric(res[i, 15:23])

        if (any(is.na(par_ss_i)) || any(is.na(par_pp_i))) return(c(ss = NA, pp = NA))

        # Reconstruct holdout split analogous to process_individual
        jag_i <- jaguar$new(id)
        track <- jag_i$get_track()
        track_cells <- jag_i$get_track_cells()

        dt_scaled <- track$dt[2:length(track$dt)] / median(na.exclude(track$dt))
        dt_discrete <- round(dt_scaled)
        outliers <- which(dt_discrete != 1)

        # Take the complement of the training set
        hold <- seq_len(ceiling(nrow(track) * 0.7))
        track_cells <- track_cells[-hold]
        outliers <- outliers[outliers > max(hold)]
        if (length(outliers) > 0) outliers <- outliers - length(hold)

        sl_emp <- na.exclude(track$sl[-outliers])
        max_dist <- ceiling(1.1 * max(sl_emp) / 1000)

        # Evaluate both models on holdout data
        obj_ss <- ss_model$prepare_objects(track_cells, max_dist, rdf = brdf)
        obj_ss$outliers <- outliers
        ll_ss <- ss_model$log_likelihood(par_ss_i, obj_ss, sim = FALSE)

        inner_size_i <- res$n_jump[i] + 1
        pp_model$propagation_steps <- max(1, ceiling(8 / inner_size_i))
        obj_pp <- pp_model$prepare_objects(track_cells, max_dist, inner_size_i,
                                           rdf = brdf)
        obj_pp$outliers <- outliers
        ll_pp <- pp_model$log_likelihood(par_pp_i, obj_pp, sim = FALSE)

        c(ss = ll_ss, pp = ll_pp)
    })

    ll_holdout <- as.data.frame(t(ll_holdout))
    ll_holdout$ID <- res$ID
    names(ll_holdout) <- c("ll_ss", "ll_pp", "ID")
    if (any(ll_holdout$ll_pp == 0)) ll_holdout$ll_pp[ll_holdout$ll_pp == 0] <- NA

    saveRDS(ll_holdout, paste0("data/output/holdout_", name, ".rds"))
}