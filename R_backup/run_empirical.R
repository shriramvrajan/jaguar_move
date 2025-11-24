rm(list = setdiff(ls(), "nbhd_full"))
source("R_backup/functions.R")     # Existing functions
source("R_backup/classes.R")       # New classes

set.seed(7)

model_type <- 2    # 1 for step selection, 2 for path propagation

config <- list(   
    # Model parameters
    model_type        = model_type,
    npar              = switch(model_type, 9, 9),    # Number of parameters
    step_size         = 1,    # Minimum step size (inner neighborhood) in pixels
    # propagation_steps = NA,    # Number of steps to simulate for path propagation (set)
    individuals       = NULL,        # NULL = all individuals, or vector of IDs

    # Holdout set parameters
    holdout_set  = FALSE,      # Whether to use holdout set for evaluation (T/F)
    holdout_frac = 0.7,       # Fraction of data to use for training

    # Parallel processing parameters
    parallel = FALSE,         # Whether to use parallel processing (T/F)
    n_cores  = 3,             # Number of cores to use if parallel

    # Model fitting options
    fit_model      = TRUE,    # Whether to fit the model (T/F)
    model_calcnull = FALSE    # Whether to calculate null model likelihood (T/F)
)

batch <- empirical_batch$new(config)
results <- batch$run_all()

saveRDS(results, paste0("data/output/empirical_results_", 
                        switch(config$model_type, "ss", "pp"), "_",
                        Sys.Date(), ".rds"))
