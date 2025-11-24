rm(list = ls())
source("R_backup/functions.R")  
source("R_backup/classes.R")       

## Autocorrelation range study
batch <- simulation_batch$new()
r1_values <- c(1, 2, 5, 7) #, 10, 15, 20, 30, 50)
batch$configs <- lapply(r1_values, function(r1) {
    simulation_config$new(
        # Simulation parameters
        name            = paste0("r1_", r1),
        obs_interval    = 2,        # Observation interval in time steps
        n_steps         = 300,     # Number of steps to simulate
        n_individuals   = 10,       # Number of individuals to simulate

        # Landscape parameters
        env_size  = 500,                   # Square side length in pixels
        autocorr_range  = r1,              # Autocorrelation range
        autocorr_strength = 5,             # Autocorrelation strength
        env_response = c(-1.5, 1.5, -0.2), # Environmental response parameters

        # Model fitting parameters
        step_size       = 1,        # Minimum step size in pixels
        n_cores         = 5         # Number of cores for parallel processing
    )
})

output <- batch$run_all(parallel =  FALSE)
saveRDS(output, paste0("simulations/results_", Sys.time(), ".rds"))

results <- batch$get_results()[[1]]
summary <- batch$get_results()[[2]]
print(summary)

batch$plot_fragmentation_effect()
