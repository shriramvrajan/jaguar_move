source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

config <- empirical_config$new(model_type = 2, parallel = TRUE, n_cores = 10)
batch <- empirical_batch$new(config)
results <- batch$run_all()