rm(list = ls())
source("R/functions.R")  
source("R/classes.R")       

test_barrier    <- TRUE
test_fragment   <- FALSE

## Barrier study
if (test_barrier) {
    batch <- simulation_batch$new()
    # b_density_values <- c(0, 15, 30, 45, 60, 75)
    b_density_values <- 50
    batch$configs <- lapply(b_density_values, function(b1) {
        simulation_config$new(
            # Simulation parameters
            name            = paste0("b1_", b1),
            obs_interval    = 10,                # Observation interval in time steps
            n_steps         = 2000,              # Number of steps to simulate
            n_individuals   = 5,                # Number of individuals to simulate

            # Environmental response parameters
            # First three for env_function(), last one is log(k_exp)
            env_response = c(-1.5, 1.5, -0.2, 0),  

            # Landscape parameters
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = 15,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            b_density         = b1,

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 5         # Number of cores for parallel processing
        )
    })
    output <- batch$run_all(parallel = TRUE)
    saveRDS(output, paste0("simulations/barriersim_", Sys.time(), ".rds"))
    batch$plot_barrier()
}

## Autocorrelation range / observation interval study
if (test_fragment) {
    batch <- simulation_batch$new()
    param_grid <- expand.grid(
        r1 = c(1, 5, 10, 20, 40), obs_interval = c(0, 1, 2, 3, 5, 8, 10),
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
            env_response = c(-1.5, 1.5, -0.2, 0),  # Environmental response parameters

            # Landscape parameters
            b_density         = 0,              # No barriers
            env_size          = 400,            # Square side length in pixels
            autocorr_range    = r1,             # Autocorrelation range
            autocorr_strength = 10,             # Autocorrelation strength
            # No barriers

            # Model fitting parameters
            step_size       = 1,        # Minimum step size in pixels
            n_cores         = 10         # Number of cores for parallel processing
        )
    })

    output <- batch$run_all(parallel = TRUE)
    saveRDS(batch, paste0("simulations/r1_obs_sweep_", Sys.time(), ".rds"))
    batch$plot_2d_sweep()
    print(batch$get_results()[[2]])
}

if (FALSE) {
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

    out0 <- readRDS("simulations/r1sim_r1_1_0_500.rds")
    out1 <- readRDS("simulations/r1sim_r1_1_1_500.rds")
    out2 <- readRDS("simulations/r1sim_r1_1_2_500.rds")

    res_table <- function(out) {
        res <- do.call(rbind, lapply(out, function(r1) {
            fits <- r1$fits
            ll_ss <- unlist(lapply(fits, function(i) {
                if (length(i$step_selection) > 1) {
                    return(i$step_selection$ll)
                } else {
                    return(NA)
                }
            }))
            ll_pp <- unlist(lapply(fits, function(i) {
                if (length(i$path_propagation) > 1) {
                    return(i$path_propagation$ll)
                } else {
                    return(NA)
                }
            }))
            return(as.data.frame(cbind(r1$config$autocorr_range, ll_ss, ll_pp)))
        }))
        names(res) <- c("r1", "ll_ss", "ll_pp")
        res$r1 <- as.numeric(res$r1)
        res$ll_ss <- as.numeric(res$ll_ss)
        res$ll_pp <- as.numeric(res$ll_pp)
        row.names(res) <- NULL
        return(res)
    }

    res0 <- res_table(out0)
    res1 <- res_table(out1)
    res2 <- res_table(out2)

    plot_res <- function(res, nm = "") {
        plotpdf(nm = paste0("figs/r1plot_", nm, ".pdf"))
        y <- tapply(res$ll_ss - res$ll_pp, res$r1, function(x) median(x, na.rm = TRUE))
        x <- sort(unique(res2$r1))
        plot(x, y, pch = 19, cex = 1.3)
        abline(h = 0, lty = 2)
        dev.off()
    }

    plot_res(res0, nm = "res0")
    plot_res(res1, nm = "res1")
    plot_res(res2, nm = "res2")
}