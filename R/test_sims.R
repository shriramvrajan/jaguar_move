source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

batch <- simulation_batch$new()$autocorr_range_study()

results <- batch$run_all(parallel = TRUE, n_cores = 5)

summary <- batch$get_performance_summary()
print(summary)

batch$plot_fragmentation_effect()