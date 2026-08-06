rm(list = ls())
source("R/functions.R")     
source("R/classes.R")       
set.seed(7)                 # For reproducibility

fit_individuals <- 0
test_holdout    <- 1
parallel        <- 1   # Whether to use parallel processing

model_type <- 1    # 1 for step selection, 2 for path propagation
env_type   <- "1o" # "1o" or "2o" or "mix", env_function argument

if (fit_individuals) {
    files <- list.files("data/output", pattern = "out_\\d+\\.rds")
    # this object is currently generated in results.R, find a better place for it
    k_fit <- readRDS("data/output/k_fitted_start_values.rds")
    ss_warm <- switch(env_type, "1o" = "emp_ss_1o_m2.rds",
                                "2o" = NULL)

    config <- list(   
        # Model parameters (npar = number of parameters)
        model_type        = model_type,
        env_type          = env_type,
        npar              = switch(env_type, "1o" = 9, "2o" = 15, "mix" = 10), 
        max_multiple      = 1,

        # jag_id$jag_id = all, or vector of specific IDs
        # individuals       = jag_id$jag_id, 
        individuals       = as.numeric(gsub("\\D", "", files)), # Save existing results

        # Holdout set parameters
        holdout_set  = TRUE,     # Whether to reserve holdout set (T/F)
        holdout_frac = 0.6,       # Proportion of data to use for training

        # Parallel processing parameters
        parallel = parallel,
        n_cores  = 10,             # Number of cores to use if parallel
        mem_budget = 0.7 * 94e9,

        # Model fitting options
        fit_model      = TRUE    # Whether to fit the model (T/F)
        # model_calcnull = FALSE    # Whether to calculate null model likelihood (T/F)
    )

    batch <- empirical_batch$new(config, k_fit = k_fit, ss_warm_par = ss_warm)
    results <- batch$run_all()

    message("Saving results...")
    saveRDS(results, paste0("data/output/emp_", 
                            switch(config$model_type, "ss", "pp"), "_",
                            config$env_type, "_m",
                            config$max_multiple, "_",
                            Sys.Date(), ".rds"))

    # source("~/memplot.R") # save memory use plot, only works on vasco
}

## Holdout set evaluation ======================================================

if (test_holdout) {

    ## Must mirror the fitting run, ideally would be read from a config object
    holdout_frac <- 0.6   # config$nholdout_frac when fitting
    max_multiple <- 1    # config$max_multiple
    scale_time <- FALSE  # config$scale_time 
    fix_depth  <- TRUE   # ?

    res <- results_set$new(
        r_ss = "data/output/emp_ss_1o_m1_holdout.rds",
        r_pp = "data/output/emp_pp_1o_m1_holdout.rds",
        env_type = env_type
    )$res_table %>% as.data.frame

    npar <- switch(env_type, "1o" = 9, "2o" = 15, "mix" = 10)
    ss_cols <- paste0("ss_par", seq_len(npar))
    pp_cols <- paste0("pp_par", seq_len(npar))

    ss_model <- step_selection_model$new()
    blank <- function(id) data.frame(ID = id, nll_ss = NA_real_, nll_pp = NA_real_,
                                   n_ss = NA_integer_, n_pp = NA_integer_)
    rows <- which(!is.na(res$ss_ll) & !is.na(res$pp_ll))

    ### Change seq_len(nrow(res)) part
    ll_holdout <- lapply(rows, function(i) {
        id <- res$ID[i]
        message(paste0("Holdout analysis for individual ", id))

        par_ss_i <- as.numeric(res[i, ss_cols])
        par_pp_i <- as.numeric(res[i, pp_cols])
        njump_i <- res$pp_njump[i]
        if (anyNA(par_ss_i) || anyNA(par_pp_i) || is.na(njump_i)) 
            return(blank(id))

        # Reconstruct holdout split analogous to process_individual
        jag_i <- jaguar$new(id, max_multiple = max_multiple, scale_time = scale_time)
        if (nrow(jag_i$track) <= 100) return(blank(id))
        track_i <- jag_i$track_cells
        hold      <- seq_len(ceiling(length(track_i) * holdout_frac))
        track_i_h <- track_i[-hold]
        out_h     <- jag_i$outliers[jag_i$outliers > max(hold)] - length(hold)

        # Step selection
        obj_ss <- ss_model$prepare_objects(track_i_h, jag_i$max_dist, rdf = brdf)
        obj_ss$outliers <- as.integer(out_h)
        ll_ss <- ss_model$log_likelihood(par_ss_i, obj_ss, sim = FALSE,
                                         env_type = env_type)
        v_ss <- setdiff(seq_along(obj_ss$obs), out_h)
        n_ss <- sum(!is.na(obj_ss$obs[v_ss]))

        # Path propagation
        setup <- jag_i$pp_objects(njump_i, env_type, holdout_frac = holdout_frac,
                                  training = FALSE)
        if (is.null(setup)) return(data.frame(ID = id, nll_ss = ll_ss, 
                                           nll_pp = NA, n_ss = n_ss, n_pp = NA))
        ll_pp <- if (!fix_depth) {
            setup$model$log_likelihood(par_pp_i, setup$objects, sim = FALSE,
                                       env_type = env_type)
        } else {
            train <- jag_i$pp_objects(njump_i, env_type,
                                      holdout_frac = holdout_frac, training = TRUE)
            r_hat <- setup$model$diagnose(par_pp_i, train$objects,
                                          env_type = env_type)$best_r
            d <- setup$model$diagnose(par_pp_i, setup$objects,
                                      env_type = env_type)
            if (r_hat + 1 > length(d$ll_by_baserate)) NA_real_
              else -d$ll_by_baserate[r_hat + 1]
        }

        v_pp <- setdiff(seq_along(setup$objects$obs), setup$objects$outliers)
        n_pp <- sum(!is.na(setup$objects$obs[v_pp]))

        data.frame(ID = id, nll_ss = ll_ss, nll_pp = ll_pp,
                   n_ss = n_ss, n_pp = n_pp)
    })

    ll_holdout <- do.call(rbind, ll_holdout)
    ll_holdout$nll_ss_per_obs <- ll_holdout$nll_ss / ll_holdout$n_ss
    ll_holdout$nll_pp_per_obs <- ll_holdout$nll_pp / ll_holdout$n_pp

    saveRDS(ll_holdout, paste0("data/output/holdout_", env_type,
                               "_m", max_multiple, "_f", holdout_frac, "_",
                               Sys.Date(), ".rds"))
}
