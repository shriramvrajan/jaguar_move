rm(list = ls())
skip_data <- TRUE
source("R/functions.R")  
source("R/classes.R")       

test_barrier    <- FALSE
test_2dsweep    <- TRUE

## Barrier density study =======================================================
if (test_barrier) {
    batch <- simulation_batch$new()
    b_density_values <- seq(0, 30, 5)
    batch$configs <- lapply(b_density_values, function(b1) {
        simulation_config$new(
            # Simulation parameters
            name            = paste0("b1_", b1),
            obs_interval    = 4,                # Observation interval in time steps
            n_steps         = 1000,              # Number of steps to simulate
            n_individuals   = 30,                # Number of individuals to simulate

            # Environmental response parameters
            # First three for env_function(), last one is log(k_exp)
            env_response = c(4, -3, 0.5, 0),  

            # Landscape parameters
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = .Machine$double.eps,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            b_density         = b1,
            b_value           = 50,
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

    plot_barrier <- function(batch_results) {
            results_df <- batch_results$results
            plot_pdf(nm = "figs/btest.pdf")
            ggplot(data = results_df, 
                aes(x = prop_barrier, y = aic_diff, col = config_name)) +
                geom_point(alpha = 0.5) +
                # geom_text(aes(label = individual, x = connect + 0.02), size = 2) +
                labs(y = "Difference in AIC (SS - PP)") +
                theme_minimal()
            dev.off()
        }
    plot_barrier(batch$get_results())
    file.remove(list.files("data/output", pattern = "sim_b1", full.names = TRUE))    
}

## Autocorrelation range / observation interval study ==========================
if (test_2dsweep) {
    batch <- simulation_batch$new()
    r1_values           <- c(.Machine$double.eps, seq(5, 25, 5))
    obs_interval_values <- 1:12
    param_grid <- expand.grid(
        r1 = r1_values, obs_interval = obs_interval_values,
        stringsAsFactors = FALSE
    )
    batch$configs <- lapply(seq_len(nrow(param_grid)), function(i) {
        r1  <- param_grid$r1[i]
        o_i <- param_grid$obs_interval[i]

        simulation_config$new(
            # Simulation parameters
            name            = paste0("r1_", r1, "obs_", o_i),
            obs_interval    = o_i,           # Observation interval in time steps
            n_steps         = 2000,          # Number of steps to simulate
            n_individuals   = 30,            # Number of individuals to simulate
            env_response = c(4, -3, 0.5, 0), # Environmental response parameters

            # Landscape parameters
            b_density         = 0,              # No barriers
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = r1,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength

            # Model fitting parameters
            step_size       = 1,      # Minimum step size in pixels
            n_cores         = 4       # Number of cores for parallel processing
        )
    })
      
    done_files <- list.files("data/output", pattern = "^sim_r1.*\\.rds$")
    batch$done <- gsub("^sim_|\\.rds$", "", done_files)
    output <- batch$run_all(parallel = TRUE)
    saveRDS(batch$get_results(), 
            paste0("simulations/r1_obs_sweep_", Sys.time(), ".rds"))  

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
        return(results)
    } 
    plot_2d_sweep(batch$get_results())
    file.remove(list.files("data/output", pattern = "sim_r1", full.names = TRUE))
}

## Other analyses ==============================================================

if (FALSE) {

    ### Plotting barrier simulation results ------------------------------------

    res1 <- readRDS("simulations/barriersim_seq_0_30_5.rds")
    res1 <- res1$results
    res1 <- res1[-which(res1$ll_step == 0 | res1$ll_pp == 0), ]
    plot_pdf(nm = "figs/btest.pdf")
    ggplot(data = res1, aes(x = b_density, y = aic_diff, col = config_name)) +
        geom_point(alpha = 0.5) +
        geom_text(aes(label = individual, x = b_density + 1), size = 2) +
        labs(y = "Difference in AIC (SS - PP)") +
        theme_minimal()
    dev.off()

    ### Plotting autocorrelation range / observation interval sweep results ----

    res2 <- readRDS("simulations/r1_obs_sweep_2026-02-13 13:20:41.184713.rds")
    res2 <- res2$summary
    # res2 <- res2[!is.na(res2$ll_diff_per_obs), ]
    plot_pdf(nm = "figs/r1test2.pdf", x = 7, y = 5)
    ggplot(data = res2, aes(x = obs_interval, y = autocorr_range, 
                            fill = ll_diff_per_obs)) +
        geom_tile() +
        # geom_text(aes(label = round(aic_diff, 2)), size = 3) +
        scale_fill_viridis_c(name = "Median LL\n(SS - PP) per obs") +
        labs(x = "Observation interval", y = "Autocorrelation range") +
        theme_minimal()
    dev.off()


    plot_pdf(nm = "test.pdf", x = 8, y = 8)
    par(mfrow = c(2, 2))
    for (i in unique(res1$config_name)) {
        subset <- res1[res1$config_name == i, ]
        hist(subset$prop_barrier, main = i)
        abline(v = median(subset$prop_barrier), col = "red", lwd = 2)
    }
    dev.off()

}