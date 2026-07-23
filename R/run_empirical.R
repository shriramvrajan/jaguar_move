rm(list = ls())
source("R/functions.R")     
source("R/classes.R")       
set.seed(7)                 # For reproducibility

fit_individuals <- 1
test_holdout    <- 0
parallel        <- 0   # Whether to use parallel processing

# 1 for step selection, 2 for path propagation
model_type <- 1
env_type   <- "1o" # "1o" or "2o" or "mix", env_function argument

if (fit_individuals) {
    files <- list.files("data/output", pattern = "out_\\d+\\.rds")
    # this object is currently generated in results.R, find a better place for it
    k_fit <- readRDS("data/output/k_fitted_start_values.rds")
    ss_warm <- switch(env_type, "1o" = readRDS("data/output/emp_ss_1o.rds"),
                                "2o" = NULL)

    config <- list(   
        # Model parameters (npar = number of parameters)
        model_type        = model_type,
        env_type          = env_type,
        npar              = switch(env_type, "1o" = 9, "2o" = 15, "mix" = 10), 

        # NULL = all individuals, or vector of specific IDs
        individuals       = jag_id$jag_id[1:3], 
        # individuals       = as.numeric(gsub("\\D", "", files)), # Save existing results

        # Holdout set parameters
        holdout_set  = TRUE,     # Whether to reserve holdout set (T/F)
        holdout_frac = 0.6,       # Proportion of data to use for training

        # Parallel processing parameters
        parallel = parallel,
        n_cores  = 20,             # Number of cores to use if parallel
        mem_budget = 0.75 * 94e9,

        # Model fitting options
        fit_model      = TRUE    # Whether to fit the model (T/F)
        # model_calcnull = FALSE    # Whether to calculate null model likelihood (T/F)
    )

    batch <- empirical_batch$new(config, k_fit = k_fit, ss_warm_par = ss_warm)
    results <- batch$run_all()

    message("Saving results...")
    saveRDS(results, paste0("data/output/emp_", 
                            switch(config$model_type, "ss", "pp"), "_",
                            Sys.Date(), ".rds"))

    source("~/memplot.R") # save memory use plot, only works on vasco
}

## Holdout set evaluation ======================================================

if (test_holdout) {
    ## NEEDS FIXING FOR MULTIPLIERS BEFORE BEING RUN AGAIN
    res <- results_set$new(
        r_ss = "data/output/empirical_ss_qcbs_holdout.rds",
        r_pp = "data/output/empirical_pp_qcbs_holdout.rds",
        env_type = env_type
    )$res_table %>% as.data.frame

    npar <- switch(env_type, "1o" = 9, "2o" = 15, "mix" = 10)
    ss_cols <- paste0("ss_par", seq_len(npar))
    pp_cols <- paste0("pp_par", seq_len(npar))

    ss_model <- step_selection_model$new()
    pp_model <- path_propagation_model$new()

    ### Change seq_len(nrow(res)) part
    ll_holdout <- sapply(seq_len(nrow(res)), function(i) {
        id <- res$ID[i]
        message(paste0("Holdout analysis for individual ", id))

        par_ss_i <- as.numeric(res[i, 3:11])
        par_pp_i <- as.numeric(res[i, 15:23])
        njump_i <- 0 # duct tape
        if (anyNA(par_ss_i) || anyNA(par_pp_i) || is.na(njump_i)) 
            return(c(ss = NA, pp = NA))

        # Reconstruct holdout split analogous to process_individual
        jag_i <- jaguar$new(id)
        track <- jag_i$get_track()
        track_cells <- jag_i$track_cells
        if (nrow(track) <= 100) return(c(ss = NA, pp = NA))

        dt_discrete <- round(track$dt[-1] / median(na.exclude(track$dt)))
        outliers_full <- which(dt_discrete != 1)
        sl_emp <- if (length(outliers_full)) track$sl[-outliers_full] else track$sl
        max_dist <- ceiling(1.1 * max(na.exclude(sl_emp)) / 1000)

        # Take the complement of the training set
        hold <- seq_len(ceiling(nrow(track) * 0.7))
        track_cells <- track_cells[-hold]
        outliers <- outliers_full[outliers_full > max(hold)] - length(hold)

        # Evaluate both models on holdout data
        obj_ss <- ss_model$prepare_objects(track_cells, max_dist, rdf = brdf)
        obj_ss$outliers <- outliers
        ll_ss <- ss_model$log_likelihood(par_ss_i, obj_ss, sim = FALSE,
                                         env_type = env_type)

        step_size_i <- njump_i + 1
        pp_model$propagation_steps <- min(max(1, ceiling(max_dist / step_size_i)), 8)
        obj_pp <- pp_model$prepare_objects(track_cells, max_dist, step_size_i,
                                           rdf = brdf)
        obj_pp$outliers <- outliers
        ll_pp <- pp_model$log_likelihood(par_pp_i, obj_pp, sim = FALSE,
                                          env_type = env_type)

        c(ss = ll_ss, pp = ll_pp)
    })

    ll_holdout <- as.data.frame(t(ll_holdout))
    names(ll_holdout) <- c("nll_ss", "nll_pp")
    ll_holdout$ID <- res$ID
    ll_holdout$nll_pp[which(ll_holdout$nll_pp == 0)] <- NA

    saveRDS(ll_holdout, paste0("data/output/holdout_", env_type, "_",
                               Sys.Date(), ".rds"))
}
