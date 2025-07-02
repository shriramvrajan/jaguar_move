source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

batch <- simulation_batch$new()$autocorr_range_study(
    # r1_values = c(1, 10, 30, 50, 80, 100),
    # n_individuals = 15,
    r1_values = 1,
    n_individuals = 1
)

for (i in seq_along(batch$configs)) {
  batch$configs[[i]]$obs_interval <- 1
  batch$configs[[i]]$n_steps <- 1000
}

results <- batch$run_all(parallel = FALSE, n_cores = 5)

# Extract performance summary
summary <- batch$get_performance_summary()
print(summary)

# Create plots
# batch$plot_fragmentation_effect()