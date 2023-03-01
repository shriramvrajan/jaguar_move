## Main script to run others

source("R/00_basics.R")

generate_data   <- FALSE
run_model       <- TRUE
analyse_output  <- FALSE

if (generate_data) {
    source("R/01_generate_data.R")
}

if (run_model) {
    source("R/02_movement_model.R")
}

if (analyse_output) {
    source("R/03_output_analysis.R")
}