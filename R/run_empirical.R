rm(list = ls())
source("R/functions.R")  # Existing functions
source("R/classes.R")       # New classes

config <- empirical_config$new(model_type = 2, parallel = FALSE, n_cores = 10)
batch <- empirical_batch$new(config)
results <- batch$run_all()
