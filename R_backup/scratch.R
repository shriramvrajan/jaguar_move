source("R/functions.R")
source("R/classes.R")

j114 <- empirical_config$new(model_type = 1, individuals = c(114))
j114_batch <- empirical_batch$new(j114)
j114_results <- j114_batch$run_all()


old <- readRDS("data/output/oldresults_empirical.rds")


      objects <- self$prepare_objects(init_point, max_dist, step_size, rdf, sim,
                                      brazil_ras)
      current <- self$build_kernel(par, objects, sim)

      # Return final probability distribution for the starting point
      return(list(
        probs = current[, 1, self$propagation_steps],  # Final step probabilities
        nbhd = objects$nbhd_i,
        par = par,
        propagation_steps = self$propagation_steps
      ))                                   