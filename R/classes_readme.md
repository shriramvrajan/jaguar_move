# Classes

## 1. path_propagation_model

propagation_steps: should be set externally

### prepare_objects

input: trajectory, rdf, sim, max_dist, step_size
output: env_i, nbhd_i, to_dest, dest, obs, outliers, mu_env, sd_env

### build_kernel()

input: par, objects, sim
output: current

### dispersal_from

input: init_point, par, step_size, n_steps, rdf, raster, threshold
output: probs, nbhd, par

### log_likelihood

input: par, objects, sim, debug
output: out (negative log-likelihood)

### fit

input: trajectory, max_dist, step_size, rdf, par_start, sim, env_raster, outliers
output: par, ll, convergence

## 2. movement_simulator

### initialize

input: config
output: movement_simulator object

### simulate_paths

input: none (based on self$config)
output: simulated paths

### fit_models

input: paths, parallel, n_cores
output: fitted models
