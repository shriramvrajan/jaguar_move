source("R/00_functions.R")  # Existing functions
source("R/classes.R")       # New classes

batch <- simulation_batch$new()$autocorr_range_study()

batch$configs <- sapply(batch$configs, function(conf) { 
    conf$obs_interval <- 5
    conf$n_steps      <- 2000
    conf$n_individuals <- 30
    return(conf)
})

results <- batch$run_all(parallel = TRUE)

summary <- batch$get_performance_summary()
print(summary)

batch$plot_fragmentation_effect()
