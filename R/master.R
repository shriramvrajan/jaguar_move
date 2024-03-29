## Main script for running others

rm(list = ls())

source("R/00_functions.R")

## Global switches =============================================================

generate_data   <- FALSE
run_model       <- TRUE
analyse_output  <- FALSE


## Set up parallel processing ==================================================

parallel        <- TRUE
ncore           <- 8
if (parallel) {
    library(doParallel)
    library(foreach)
    registerDoParallel(ncore)
    message(paste0("number of workers: ", getDoParWorkers()))
}

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
      refit_model0  <- FALSE            # Refit traditional SSF model
    holdout_set     <- TRUE             # Hold out a set of locations
    model_calcnull  <- FALSE            # Calculate null likelihoods 
                                        # refit_model must be TRUE for this one
                                    
    ## Parameters                                    
    npar            <- 7              # Number of parameters in current model
    sim_steps       <- 10             # How many steps to simulate forward
    buffersize      <- 1              # Jaguars move 1px (1km) at a time
    n_iter          <- nrow(jag_id)   # Number of individuals
    holdout_frac    <- 0.7            # Fraction of points to use for fitting
    
    message("Parameters set")
    message(paste0("Number of jaguars: ", n_iter))
    message(paste0("Holdout fraction: ", holdout_frac))
    message(paste0("Number of parameters: ", npar))
    message(paste0("Number of simulation steps: ", sim_steps))
    message(paste0("Buffer size: ", buffersize))
    message("=========================================")

    outfiles <- list.files("data/output")
    done <- gsub(".RDS", "", outfiles) %>% 
        gsub("par_out_", "", .) %>% 
        gsub("NA_", "", .) %>%
        as.numeric()
    done <- done[!is.na(done)]
    i_todo <- setdiff(seq_len(n_iter), done)
    message(paste0("Number of jaguars to process: ", length(i_todo)))

    source("R/02_movement_model.R")
}

## Output analysis =============================================================

if (analyse_output) {
    source("R/03_output_analysis.R")
}
