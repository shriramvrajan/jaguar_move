## Main script to run others

source("R/00_basics.R")

generate_data   <- FALSE
run_model       <- TRUE
analyse_output  <- FALSE

if (generate_data) {
    source("R/01_generate_data.R")
}

if (run_model) {
    i_initial       <- 1              # Individual to start at
    buffersize      <- 1              # Jaguars move 1px (1km) at a time
    model_calcnull  <- TRUE           # Calculate null likelihoods?
    refit_homes     <- FALSE          # Refit home ranges (AKDE) 
    refit_turns     <- FALSE          # Refit turn distributions (MM)
    # model_fiteach   <- FALSE
    source("R/03_movement_model.R")
}

if (analyse_output) {
    source("R/04_output_analysis.R")
}