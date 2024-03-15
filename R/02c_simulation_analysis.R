source("R/00_functions.R")

## Switches ====================================================================
plot_aic      <- F
gen_parscape  <- T
plot_parscape <- F

parallel_setup(20)

## Load data ===================================================================
message("Loading data")
params <- readRDS("data/output/simulations/params.rds")
sim_interval <- params[[4]]
sim_n <- params[[7]]
paths <- readRDS("data/output/simulations/paths.rds")
env01 <- rast("data/output/simulations/env01.tif")
env01 <- list(env01, raster_to_df(env01))
fit <- load_if_exists(paste0("par_out_", 1:sim_n, ".rds"), 
                             dir = "data/output/simulations") %>%
        do.call(rbind, .) %>% 
        as.data.frame() 
fit$id <- seq_len(nrow(fit))
names(fit) <- c("m1", "c1", "b1", "a1", "m2", "c2", "b2", "a2", "id")
posna <- which(is.na(fit$m1)) # positions of NA
fit <- fit[-posna, ]

## Plot model AIC ==============================================================
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
    aic1 <- 2 * 4 - 2 * ll1
    aic2 <- 2 * 4 - 2 * ll2
    # CHECK NUMBER OF PARAMETERS ^
    par(mfrow = c(1, 1))
    plot(aic1, aic2)
    abline(0, 1)
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

# Generate optim landscape =====================================================
if (gen_parscape) {
    message("Generating fitting landscape...")

    x <- seq(-4, 4, length.out = 50)
    y <- seq(-4, 4, length.out = 50)
    testvals <- expand.grid(x, y)
    names(testvals) <- c("x", "y")

    to_gen <- seq_len(nrow(fit))[-which(is.na(fit$a1))]
    for (i in to_gen) {
    # foreach(i = to_gen) %dopar% {    
        message(paste0("Generating landscape", i, " / ", nrow(fit)))
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

        message("Calculating log likelihoods")
        testvals$ll <- unlist(foreach(j = seq_len(nrow(testvals)), .export = c("dest")) %dopar% {
            message(paste0("Fitting row #: ", j, " / ", nrow(testvals)))
            val <- as.numeric(testvals[j, ])
            ll <- log_likelihood(val, objects1)
            message(paste0("Fitted row #: ", j, " / ", nrow(testvals)))
            return(ll)
        })
        message("Aggregating results")
        saveRDS(testvals, "data/output/simulations/testll.rds")
    }
}

## Plot parameter landscape ====================================================
if (plot_parscape) {
        testvals <- readRDS("data/output/simulations/testll.rds")
        plot <- ggplot(testvals, aes(x = x, y = y, fill = ll)) +
            geom_tile() +
            scale_fill_viridis_c(option = "turbo") +
            theme_minimal() +
            labs(title = paste0("Log likelihood surface", i),
                x = "a1", y = "b1") +
            theme(legend.position = "none") +
            geom_point(aes(x = 3, y = -2), color = "magenta", size = 3)
        ggsave(paste0("data/output/simulations/ll_surface_", i, ".png"), plot,
               bg = "white")
}


x1 <- seq(0, 7, length.out = 100)

y0 <- 1 / (1 + exp(2 - 0.2 * x1 - 0.2 * x1^2)) # real parameter values, 0 -1 -1
y1 <- lapply(fit$id, function(i) {
    out <- 1 / (1 + exp(fit$c1[i] + fit$b1[i] * x1 + fit$a1[i] * x1^2))
})

points <- lapply(fit$id, function(i) {
    path <- paths[[i]]
    path$move <- c(0, sqrt(diff(path$x)^2 + diff(path$y)^2))
    path$move[path$move > 0] <- 1
    path$env <- env01[[2]]$sim1[cellFromXY(env01[[1]], path[, 1:2])]
    return(path[, c("move", "env")])
})

par(mfrow = c(4, 5))
plot(x1, y0, type = "l", main = "Generating function", ylim = c(0, 1))
for (i in seq_len(nrow(fit))) {
    plot(x1, y1[[i]], col = "#363636", type = "l", ylim = c(0, 1),
         main = paste0("Fitted individual ", i))
    points(points[[i]]$env, points[[i]]$move, col = rgb(0, 0, 0, 0.3), 
           cex = 0.5, pch = 19)
}

