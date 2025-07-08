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

# Manual debugging - extract one trajectory and compare objects
config <- batch$configs[[1]]
landscape <- batch$results[[names(batch$results)[1]]]$landscape
paths <- batch$results[[names(batch$results)[1]]]$paths

# Prepare trajectory (copy from fit_models)
jag_traject <- lapply(paths, function(p) {
    out <- cbind(p$x, p$y)
    ind <- seq(1, nrow(out), config$obs_interval + 1)
    out <- out[ind, ]
    return(out)
})

jag_traject_cells <- lapply(jag_traject, function(tr) {
    out <- terra::cellFromRowCol(landscape$raster, tr[, 1], tr[, 2])
    return(out)
})

rdf <- raster_to_df(landscape$raster)
env_ras <- landscape$raster
nbhd_full <<- make_nbhd(i = seq_len(ncell(landscape$raster)), 
                    sz = config$step_size, r = landscape$raster, rdf = rdf)

# Manual object preparation
trajectory <- jag_traject_cells[[1]]  # First individual
max_dist <- config$max_dist()
sim_steps <- config$sim_steps()

ss_model <- step_selection_model$new()
pp_model <- path_propagation_model$new()

par_start <- c(-1.5, 1.5, -0.2)

# Fit both models
ss_result <- ss_model$fit(trajectory, max_dist, rdf, par_start, 
                        sim = TRUE, env_raster = landscape$raster)
pp_result <- pp_model$fit(trajectory, max_dist, rdf, par_start,
                        sim = TRUE, env_raster = landscape$raster, 
                        sim_steps = sim_steps)

# Compare objects side by side
ss_objects <- ss_model$prepare_objects(trajectory, max_dist, rdf, sim = TRUE, env_raster = landscape$raster)
pp_objects <- pp_model$prepare_objects(trajectory, max_dist, rdf, sim = TRUE, env_raster = landscape$raster, sim_steps = sim_steps)