rm(list = ls())
skip_data <- TRUE
source("R/functions.R")  
source("R/classes.R")       

test_barrier    <- TRUE
test_2dsweep    <- FALSE

## Plotting functions ==========================================================

plot_2d_sweep <- function(batch_results) {
    results <- batch_results$summary

    # 2D heatmap of median LL difference per observation
    p <- ggplot(results, aes(x = obs_interval, y = autocorr_range,
                                fill = ll_diff_per_obs)) +
        geom_tile() +
        geom_text(aes(label = round(ll_diff_per_obs, 3)), size = 3) +
        scale_fill_viridis_c(name = "Median LL\n(SS - PP) per obs") +
        labs(x = "Observation interval", y = "Autocorrelation range") +
        theme_minimal()
    ggsave(paste0("figs/sims/2d_sweep_", Sys.time(), ".pdf"), p,
            width = 8, height = 6)
    print(p)
    return(summary_df)
}

plot_barrier <- function(batch_results) {
        results_df <- batch_results$results
        # AIC difference vs barrier density
        lld <- tapply(results_df$ll_diff_per_obs, results_df$b_density, 
                        function(x) median(x, na.rm = T))
        plot_pdf(nm = paste0("figs/sims/barrier_study", Sys.time(), ".pdf"), 
            x = 6, y = 4)
        plot(results_df$connect, results_df$aic_diff,
            xlab = "Landscape connectivity", 
            ylab = "Difference in AIC (SS - PP)",
            pch = 19, cex = 1, col = rgb(0, 0, 0, 0.5))
        dev.off()          
        return(results_df)
}

## Barrier density study =======================================================
if (test_barrier) {
    batch <- simulation_batch$new()
    b_density_values <- 50 # seq(10, 100, 5)
    batch$configs <- lapply(b_density_values, function(b1) {
        simulation_config$new(
            # Simulation parameters
            name            = paste0("b1_", b1),
            obs_interval    = 3,                # Observation interval in time steps
            n_steps         = 1000,              # Number of steps to simulate
            n_individuals   = 1,                # Number of individuals to simulate

            # Environmental response parameters
            # First three for env_function(), last one is log(k_exp)
            env_response = c(4, -3, 0.5, 0),  

            # Landscape parameters
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = 5,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            b_density         = b1,
            b_value           = 99,
            b_width           = 2,

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 6         # Number of cores for parallel processing
        )
    })
          
    done_files <- list.files("data/output", pattern = "^sim_b1.*\\.rds$")
    batch$done <- gsub("^sim_|\\.rds$", "", done_files)
    output <- batch$run_all(parallel = TRUE)
    saveRDS(batch$get_results(), 
        paste0("simulations/barriersim_", Sys.time(), ".rds"))
    plot_barrier(batch$get_results())
    file.remove(list.files("data/output", pattern = "sim_b1", full.names = TRUE))    
}

## Autocorrelation range / observation interval study ==========================
if (test_2dsweep) {
    batch <- simulation_batch$new()
    param_grid <- expand.grid(
        r1 = 1:100, obs_interval = 1:20,
        stringsAsFactors = FALSE
    )
    batch$configs <- lapply(seq_len(nrow(param_grid)), function(i) {
        r1  <- param_grid$r1[i]
        o_i <- param_grid$obs_interval[i]

        simulation_config$new(
            # Simulation parameters
            name            = paste0("r1_", r1, "obs_", o_i),
            obs_interval    = o_i,                # Observation interval in time steps
            n_steps         = 2000,              # Number of steps to simulate
            n_individuals   = 30,                # Number of individuals to simulate
            env_response = c(4, -3, 0.5, 0),     # Environmental response parameters

            # Landscape parameters
            b_density         = 0,              # No barriers
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = r1,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 6         # Number of cores for parallel processing
        )
    })
      
    done_files <- list.files("data/output", pattern = "^sim_r1.*\\.rds$")
    batch$done <- gsub("^sim_|\\.rds$", "", done_files)
    output <- batch$run_all(parallel = TRUE)
    saveRDS(batch$get_results(), 
            paste0("simulations/r1_obs_sweep_", Sys.time(), ".rds"))   
    plot_2d_sweep(batch$get_results())
    # file.remove(list.files("data/output", pattern = "sim_r1", full.names = TRUE))
}
