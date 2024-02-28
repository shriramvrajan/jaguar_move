rerun_sim <- FALSE
ifelse(rerun_sim, source("R/01b_simulations.R"), source("R/00_functions.R"))

plot_aic      <- FALSE
gen_parscape  <- TRUE
plot_parscape <- FALSE

message("Loading data")
params <- readRDS("data/output/simulations/params.RDS")
sim_interval <- params[[5]]
sim_n <- params[[8]]

paths <- readRDS("data/output/simulations/paths.RDS")

env01 <- rast("data/output/simulations/env01.tif")
env01 <- list(env01, raster_to_df(env01))

fit <- load_if_exists(paste0("par_out_", 1:sim_n, ".RDS"), 
                             dir = "data/output/simulations") %>%
        do.call(rbind, .) %>% 
        as.data.frame() 
names(fit) <- c("a1", "b1", "a2", "b2")

if (plot_aic) {
    par(mfrow = c(1, 2))
    plot(fit$a1, fit$b1)
    abline(h = 0)
    abline(v = 0)
    plot(fit$a2, fit$b2)
    abline(h = 0)
    abline(v = 0)

    sim_n <- nrow(fit)
    ll1 <- -unlist(load_if_exists(paste0("ll_fit1", 1:sim_n, ".rds"), 
                                dir = "data/output/simulations"))
    ll2 <- -unlist(load_if_exists(paste0("ll_fit2", 1:sim_n, ".rds"), 
                                dir = "data/output/simulations"))
    aic1 <- 2 * 2 - 2 * ll1
    aic2 <- 2 * 2 - 2 * ll2

    par(mfrow = c(1, 2))
    plot(aic1, aic2)
    abline(0, 1)
    hist(aic1 - aic2, breaks = 20)
    abline(v = 0, lty = 2)
}

# Prepare simulation data ======================================================
message("Preparing simulation data...")
envdf <- env01[[2]]
env_index <- seq_len(nrow(envdf))
envdf <- cbind(envdf, env_index)
jag_traject <- lapply(paths, function(p) {
    out <- cbind(p$x, p$y)
    ind <- seq(1, nrow(out), sim_interval)
    out <- out[ind, ]
    return(out)
})
jag_traject_cells <- lapply(jag_traject, function(tr) {
    raster::cellFromXY(env01[[1]], tr[, 1:2])
})
dist <- lapply(jag_traject, function(tr) {
    out <- c(0, sqrt(diff(tr[, 1])^2 + diff(tr[, 2])^2))
    return(out)
})
max_dist <- ceiling(max(unlist(dist)) * 1.5)
step_range <- (2 * max_dist + 1) ^ 2
nbhd0 <- make_nbhd(i = seq_len(nrow(env01[[2]])), sz = buffersize, 
                   r = env01[[1]], rdf = env01[[2]]) 

# Generating optim landscape ===================================================
if (gen_parscape) {
    message("Generating fitting landscape...")

    x <- seq(-4, 4, length.out = 50)
    y <- seq(-4, 1, length.out = 50)
    testvals <- expand.grid(x, y)
    names(testvals) <- c("x", "y")

    i <- 20
    traject <- jag_traject_cells[[i]]
    prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
                    rdf = env01[[2]])
    env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
        # Normalizing desired environmental variables for extended neighborhood
    env1 <- env1[nbhd_index]
    # Make indexing consistent with env
    names(env1) <- seq_len(length(nbhd_index))
    sim_steps <- sim_interval * 2

    objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
    parallel_setup(20)

    message("Calculating log likelihoods")
    testvals$ll <- unlist(foreach(j = seq_len(nrow(testvals)), .export = c("dest")) %dopar% {
        message(paste0("Fitting row #: ", j, " / ", nrow(testvals)))
        val <- as.numeric(testvals[j, ])
        ll <- log_likelihood(val, objects1)
        message(paste0("Fitted row #: ", j, " / ", nrow(testvals)))
        return(ll)
    })
    message("Aggregating results")
    saveRDS(testvals, "data/output/simulations/testll.RDS")
}

if (plot_parscape) {
    testvals <- readRDS("data/output/simulations/testll.RDS")
    plot <- ggplot(testvals, aes(x = x, y = y, fill = ll)) +
        geom_tile() +
        scale_fill_viridis_c() +
        theme_minimal() +
        labs(title = "Log likelihood surface",
            x = "a1", y = "b1") +
        theme(legend.position = "none")
    plot
}