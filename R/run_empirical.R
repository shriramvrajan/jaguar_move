rm(list = ls())
source("R/functions.R")     # Existing functions
source("R/classes.R")       # New classes

set.seed(1)
model_type <- 2  # 1: step selection, 2: path propagation
config <- empirical_config$new(model_type = model_type,
                               parallel = TRUE, n_cores = 6)

batch <- empirical_batch$new(config)

results <- batch$run_all()
saveRDS(results, paste0("data/output/empirical_results_", 
                        switch(model_type, "ss", "pp"), "_",
                        Sys.Date(), ".rds"))
