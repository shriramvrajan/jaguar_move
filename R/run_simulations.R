rm(list = ls())
skip_data <- TRUE
source("R/functions.R")  
source("R/classes.R")       

test_jumpgrain  <- 1
test_barrier    <- 0
test_2dsweep    <- 0

## Jump x grain test ===========================================================

# Hypothesis: ll_diff per obs depends on the landscape through 
#   p = obs_interval + 1            i.e. substeps per fix
#   s = gen_step / autocorr_range   i.e. substep in units of grain (?)
if (test_jumpgrain) {
    param_grid <- expand.grid(gen_step = c(1, 2, 4),
                              s        = c(0.05, 0.1, 0.2, 0.4),
                              p        = c(2, 4, 8))
    param_grid$autocorr_range <- param_grid$gen_step / param_grid$s
    param_grid$obs_interval   <- param_grid$p - 1
    param_grid <- param_grid[param_grid$autocorr_range <= 40, ]  # keep patches < env_size/10

    batch <- simulation_batch$new()
    batch$configs <- lapply(seq_len(nrow(param_grid)), function(i) {
        g <- param_grid[i, ]
        simulation_config$new(
            name           = sprintf("jg_j%s_s%s_p%s", g$gen_step, g$s, g$p),
            obs_interval   = g$obs_interval,
            gen_step       = g$gen_step,   # truth
            step_size      = g$gen_step,   # model d_jump: well-specified
            autocorr_range = g$autocorr_range,
            n_steps = 1000, n_individuals = 12,
            env_response = c(4, -3, 0.5, 0),
            b_density = 0, env_size = 400, autocorr_strength = 10,
            n_cores = 4)
    })

    done_files <- list.files("data/output", pattern = "^sim_jg.*\\.rds$")
    batch$done <- gsub("^sim_|\\.rds$", "", done_files)
    batch$run_all(parallel = TRUE)
    out <- batch$get_results()
    saveRDS(out, paste0("simulations/jumpgrain_",
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))

    r    <- out$summary
    r$p  <- r$obs_interval + 1
    r$s  <- r$gen_step / r$autocorr_range

    pc <- ggplot(r, aes(s, ll_diff_per_obs, colour = factor(p))) +
        geom_line(aes(group = interaction(p, gen_step)), alpha = 0.35) +
        geom_point(aes(shape = factor(gen_step)), size = 2.5) +
        scale_x_log10() +
        labs(x = "substep length / autocorrelation range  (s)",
             y = "median LL (SS - PP) per obs",
             colour = "substeps per fix (p)", shape = "jump (cells/tick)") +
        theme_minimal()
    ggsave(paste0("figs/sims/jumpgrain_collapse_",
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"), pc, width = 8, height = 5)
    print(pc)
    file.remove(list.files("data/output", pattern = "sim_jg", full.names = TRUE))
}

## Barrier density study =======================================================
if (test_barrier) {
    batch            <- simulation_batch$new()
    b_density_values <- seq(0, 30, 5)
    batch$configs    <- lapply(b_density_values, function(b1) {
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
            autocorr_range    = 1,             # Autocorrelation range
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
        paste0("simulations/barriersim_", 
            format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))    

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
    r1_values           <- c(0.01, seq(5, 35, 5))
    obs_interval_values <- seq(0:10)
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
            n_steps         = 1000,          # Number of steps to simulate
            n_individuals   = 12,            # Number of individuals to simulate
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
            paste0("simulations/r1_obs_sweep_", 
                format(Sys.time(), "%Y%m%d_%H%M%S"), ".rds"))  

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
        ggsave(paste0("figs/sims/2d_sweep_",
                format(Sys.time(), "%Y%m%d_%H%M%S"), ".pdf"), p,
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

## scratch

ww <- readRDS("simulations/r1_obs_sweep_20260713_143155.rds")$results
plot(ww$autocorr_range, ww$ll_diff_per_obs)
abline(h = 0)
plot(tapply(ww$ll_diff_per_obs, ww$obs_interval, mean))
