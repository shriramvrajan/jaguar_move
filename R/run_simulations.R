rm(list = ls())
source("R_backup/functions.R")  
source("R_backup/classes.R")       

test_barrier    <- FALSE
test_fragment   <- TRUE

## Barrier study
if (test_barrier) {
    batch <- simulation_batch$new()
    b_density_values <- c(0, 10, 20, 30, 40, 50)
    batch$configs <- lapply(b_density_values, function(b) {
        simulation_config$new(
            # Simulation parameters
            name            = paste0("b_", b),
            obs_interval    = 2,                # Observation interval in time steps
            n_steps         = 2000,              # Number of steps to simulate
            n_individuals   = 10,                # Number of individuals to simulate
            env_response = c(-1.5, 1.5, -0.2),  # Environmental response parameters

            # Landscape parameters
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = 15,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            b_density         = b,

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 5         # Number of cores for parallel processing
        )
    })
    output <- batch$run_all(parallel = FALSE)
    saveRDS(output, paste0("simulations/barriersim_", Sys.time(), ".rds"))

    results <- batch$get_results()[[1]]
    summary <- batch$get_results()[[2]]
    print(summary)

}

## Autocorrelation range / habitat fragmentation study
if (test_fragment) {
    batch <- simulation_batch$new()
    # r1_values <- c(1, 3, 5, 10, 25, 20, 25, 30, 40, 50)
    r1_values <- 1
    batch$configs <- lapply(r1_values, function(r1) {
        simulation_config$new(
            # Simulation parameters
            name            = paste0("r1_", r1),
            obs_interval    = 0,                # Observation interval in time steps
            n_steps         = 500,              # Number of steps to simulate
            n_individuals   = 10,                # Number of individuals to simulate
            env_response = c(-1.5, 1.5, -0.2),  # Environmental response parameters

            # Landscape parameters
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = r1,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            # No barriers

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 7         # Number of cores for parallel processing
        )
    })

    output <- batch$run_all(parallel = FALSE)
    saveRDS(output, paste0("simulations/r1sim_", batch$configs[[1]]$name, "_", 
            batch$configs[[1]]$obs_interval, "_", batch$configs[[1]]$n_steps, ".rds"))

    results <- batch$results[[1]]$fits
    diff <- sapply(seq_len(batch$configs[[1]]$n_individuals), function(i) {
        print(i)
        fit <- results[[i]]
        if (length(fit$path_propagation) == 1) {
            return(NA)
            next
        }
        return(fit$path_propagation$ll - fit$step_selection$ll)
    })

    plotpdf(nm = paste0("figs/hist", batch$configs[[1]]$name, ".pdf"))
    hist(diff, 10, col = rgb(0, 0, 1, 0.3), border = NA, main = paste(length(which(diff < 0)), length(which(diff >0))))
    abline(v = 0, lty = 2)
    dev.off()
}

plot_parameter_recovery <- function(fits) {
    par_ss <- lapply(fits, function(i) i$step_selection$par) 
    par_pp <- lapply(fits, function(i) i$path_propagation$par)
    par0 <- c(-1.5, 1.5, -0.2)
    x <- seq(0, 10, 0.1)
    par(mfrow = c(1, 2))
    plot_curve(par0)
    for (i in seq_len(length(par_ss))) {
        lines(x, plot_curve(par_ss[[i]][1:3], values = TRUE), col = "blue")
    }
    plot_curve(par0)
    for (i in seq_len(length(par_pp))) {
        lines(x, plot_curve(par_pp[[i]][1:3], values = TRUE), col = "blue")
    }
}

bat <- readRDS("simulations/r1sim_r1_1_1_1000.rds")
plotpdf(x = 8, y = 4)
plot_parameter_recovery(bat$r1_1$fits)
dev.off()
