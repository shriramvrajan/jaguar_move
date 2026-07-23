source("R/functions.R")
source("R/classes.R")

id  <- 85

jag <- jaguar$new(85)

ss_res <- readRDS("data/output/emp_ss_1o.rds")
par0   <- ss_res[[which(jag_id$jag_id == id)]]$par
result <- jag$fit_pp(par0, env_type = "1o")
