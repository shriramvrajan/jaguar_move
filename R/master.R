## Main script for running others
if (!exists("loaded")) {
    source("R/00_functions.R")
} else {
    message("Functions already loaded")
}
loaded <- TRUE

## Global switches =============================================================

generate_data   <- FALSE
run_model       <- TRUE
analyse_output  <- FALSE

debug_02        <- FALSE

set.seed(12345)

## Data generation =============================================================

if (generate_data) {
    source("R/01_generate_data.R")
}

## Movement model ==============================================================

if (run_model) {
    
    ## Not important right now 
    refit_homes     <- FALSE            # Refit home ranges (AKDE) 
    refit_turns     <- FALSE            # Refit turn distributions (MM)

    ## Actual fitting options
    refit_model     <- TRUE             # Refit movement model parameters
    model_type      <- 2                # 1 = tradSSF, 2 = path propagation
    holdout_set     <- FALSE             # Hold out a set of points
    # holdout_frac    <- 0.7              # Fraction of points to use for fitting
    model_nofit     <- FALSE            # Don't fit the model, just simulate
        model_calcnull  <- FALSE            # Calculate null likelihoods 
                                        # refit_model must be TRUE for this one
    ## Parameters                                    
    npar            <- 7              # Number of parameters in current model 
    step_size       <- 1              # Jaguars move up to n px (n km) at a time
    sim_steps       <- 6              # How many steps to simulate forward
    i_override      <- NULL           
    # Which jaguars, set to NULL to fit all
    
    ## Set up parallel processing
    parallel        <- TRUE
    ncore           <- switch(model_type, 12, 7)
    if (parallel) {
        library(doParallel)
        library(foreach)
        registerDoParallel(ncore)
        message(paste0("number of workers: ", getDoParWorkers()))
    }

    message("Parameters set")
    message(paste0("Number of parameters: ", npar))
    message(paste0("Number of simulation steps: ", sim_steps))
    message(paste0("Step size: ", step_size))
    if (holdout_set) message(paste0("Holdout fraction: ", holdout_frac))
    message("============================================")

    if (!is.null(i_override)) {
        i_todo <- i_override
    } else {
        outfiles <- list.files("data/output")
        done <- gsub(".rds", "", outfiles) %>% 
            gsub("out_", "", .) %>% 
            gsub("ll_null_", "", .) %>%
            gsub("ll_", "", .) %>%
            gsub("NA_", "", .) %>%
            as.numeric()
        done <- done[!is.na(done)]
        i_todo <- setdiff(seq_len(nrow(jag_id)), done)
    }
    message(paste0("Number of jaguars to process: ", length(i_todo)))

    if (!exists("nbhd0")) {
        message(paste("Generating", (step_size * 2 + 1)^2, 
                      "cell neighborhood for each cell in Brazil"))
        nbhd0 <- make_nbhd(i = seq_len(nrow(brdf)), sz = step_size) 
    } else {
        message("Global neighborhood matrix already generated")
    }

    source("R/02_movement_model.R")
}

## Output analysis =============================================================

if (analyse_output) {
    source("R/03_output_analysis.R")
}
