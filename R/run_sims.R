source("R/functions.R")  # Existing functions
source("R/classes.R")       # New classes

batch <- simulation_batch$new()$autocorr_range_study()

batch$configs <- sapply(batch$configs, function(conf) { 
    conf$obs_interval  <- 5
    conf$n_steps       <- 2000
    conf$n_individuals <- 50
    conf$n_cores       <- 5 
    return(conf)
})

output <- batch$run_all(parallel =  FALSE)
# saveRDS(results, paste0("simulations/results_", Sys.time(), ".rds"))

results <- batch$get_results()[[1]]
summary <- batch$get_results()[[2]]
print(summary)

batch$plot_fragmentation_effect()
