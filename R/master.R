## Main script to run others
rm(list = ls())

source("R/00_functions.R")

generate_data   <- FALSE
run_model       <- TRUE
analyse_output  <- FALSE

if (generate_data) {
    source("R/01_generate_data.R")
}

if (run_model) {
    
    ## Switches 

    refit_homes     <- FALSE            # Refit home ranges (AKDE) 
    refit_turns     <- FALSE            # Refit turn distributions (MM)
    refit_model     <- TRUE             # Refit movement model parameters
    model_calcnull  <- FALSE            # Calculate null likelihoods 
                                        # refit_model must be TRUE for this one
    refit_model0    <- FALSE            # Refit traditional SSF model
                                        
    npar            <- 7              # Number of parameters in current model
    sim_steps       <- 25             # How many steps to simulate forward

    i_initial       <- 1              # Individual to start at
    buffersize      <- 1              # Jaguars move 1px (1km) at a time
    n_iter          <- nrow(jag_id)   # Number of individuals

    source("R/02_movement_model.R")
}

if (analyse_output) {
    source("R/03_output_analysis.R")
}