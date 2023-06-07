source("R/00_functions.R")

paths <- readRDS("data/output/simulations/paths.RDS")
probs <- load_if_exists(paste0("p", 1:sim_n, ".RDS"), 
                        dir = "data/output/simulations")
ll <- load_if_exists(paste0("ll", 1:sim_n, ".RDS"), 
                       dir = "data/output/simulations")

env01 <- readRDS("data/output/simulations/env01.RDS")
env02 <- readRDS("data/output/simulations/env02.RDS")

res <- lapply(1:sim_n, function(i) {
    path <- paths[[i]]
    prob <- probs[[i]]
    ll  <- ll[[i]]

    state <- path$state[seq(1, nrow(path), sim_interval)]
    p     <- prob[5, ]

    hist(p[state == 1], 30, main = paste0("Path #", i), col = rgb(0, 0, 1, 0.5),
         xlim = c(0, max(p, na.rm = TRUE)), border = NA)
    hist(p[state == 2], 30, add = TRUE, col = rgb(1, 0, 0, 0.5), border = NA)

    x <- readline("Press [enter] to continue") 
})
