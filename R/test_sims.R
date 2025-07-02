source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

batch <- simulation_batch$new()$create_fragmentation_study(
    r1_values = c(1, 10, 30, 50, 80, 100),
    n_individuals = 5
)
for (i in seq_along(batch$configs)) {
  batch$configs[[i]]$obs_interval <- 10
}

results <- batch$run_all(parallel = TRUE, n_cores = 5)

# Extract performance summary
summary <- batch$get_performance_summary()
print(summary)

# Create plots
batch$plot_fragmentation_effect()

