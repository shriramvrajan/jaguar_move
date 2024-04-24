source("R/00_functions.R")

## Switches ====================================================================
plot_aic    <- F
param_plots <- T
debug_fit   <- F

simdir      <- "simulations/s18/"
parallel_setup(1)

## Load data ===================================================================
message("Loading data")
params       <- readRDS(paste0(simdir, "params.rds"))   
sim_interval <- params$sim_interval
sim_n        <- params$sim_n
paths <- readRDS(paste0(simdir, "paths.rds"))
env01 <- rast(paste0(simdir, "env01.tif"))
env01 <- list(env01, raster_to_df(env01))
print(params)

fit <- load_if_exists(paste0("par_out_", 1:sim_n, ".rds"), dir = simdir) %>%
        do.call(rbind, .) %>% 
        as.data.frame() 
fit$id <- seq_len(nrow(fit))
if (ncol(fit) == 9) {
    names(fit) <- c("m1", "c1", "b1", "a1", "m2", "c2", "b2", "a2", "id")
} else if (ncol(fit) == 7) {
    names(fit) <- c("c1", "b1", "a1", "c2", "b2", "a2", "id")
} else if (ncol(fit) == 4) {
    names(fit) <- c("c1", "b1", "a1", "id")
} else {
    stop("Unexpected number of columns in fit")
}

if (any(is.na(fit$c1))) {
    posna <- which(is.na(fit$c1)) # positions of NA
    fit <- fit[-posna, ]
}

## Prepare simulation data =====================================================
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

## Explore parameter values ====================================================
x1 <- seq(0, 8, length.out = 100)
par0 <- unlist(params[10:12])

# generating parameter values
y0 <- 1 / (1 + exp(par0[1] + par0[2] * x1 + par0[3] * x1^2)) 
# fitted parameter values
yhat <- 1 / (1 + exp(median(fit$c1) + median(fit$b1) * x1 + median(fit$a1) * x1^2))
yhat2 <- 1 / (1 + exp(mean(fit$c1) + mean(fit$b1) * x1 + mean(fit$a1) * x1^2))
y1 <- lapply(seq_len(nrow(fit)), function(i) {
    out <- 1 / (1 + exp(fit$c1[i] + fit$b1[i] * x1 + fit$a1[i] * x1^2))
})
# y2 <- lapply(seq_len(nrow(fit)), function(i) {
#     out <- 1 / (1 + exp(fit$c2[i] + fit$b2[i] * x1 + fit$a2[i] * x1^2))
# })

points <- lapply(fit$id, function(i) {
    path <- paths[[i]]
    path$move <- c(0, sqrt(diff(path$x)^2 + diff(path$y)^2))
    path$move[path$move > 0] <- 1
    path$env <- env01[[2]]$sim1[cellFromXY(env01[[1]], path[, 1:2])]
    return(path[, c("move", "env")])
})

if (param_plots) {
    # Generated + fitted, all on same plot
    par(mfrow = c(1, 1))
    plot(x1, y0, type = "l", lwd = 3, ylim = c(0, 1))
    lines(x1, yhat, lwd = 3, col = "blue")
    # lines(x1, yhat2, lwd = 3, col = rgb(1, 0, 0, 0.8))
    for (i in seq_len(nrow(fit))) {
        lines(x1, y1[[i]], col = rgb(0, 0, 1, 0.3), lwd = 1)
        # lines(x1, y2[[i]], col = rgb(1, 0, 0, 0.5), lwd = 1.5)
        env <- points[[i]]$env
        move <- points[[i]]$move
        mw <- moving_window(env, move, window = 0.5)
        # lines(mw$x, mw$y, col = rgb(0, 0, 0, 0.2), lwd = 1)
        # abline(h = exp01(params$par_move), lty = 2, col = "red")
        # model <- glm(move ~ env, family = binomial)
    }
}

mm <- unlist(lapply(paths, function(x) mean(x$att, na.rm = T)))


## Debug =======================================================================

if (debug_fit) {
    # What was I trying to do here...
    parallel_setup(ncore_fit)
    llike <- load_if_exists(paste0(simdir, "ll_fit1", 1:sim_n, ".rds"), 
                            dir = ".") %>%
             unlist(.)
    par_fitted <- load_if_exists(paste0(simdir, "par_out_", 1:sim_n, ".rds"), 
                                 dir = ".") %>%
                  do.call(rbind, .) %>%
                  as.data.frame()
    # out <- foreach(i = seq_along(llike), .combine = c) %dopar% {
    out <- sapply(seq_along(llike), function(i) {
        ## !! JUST TESTING WITH ONE FIRST, generalize later !!
        i <- 12
        traject <- jag_traject_cells[[i]]
        prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
                rdf = env01[[2]])
        env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
        # Normalizing desired environmental variables for extended neighborhood
        env1 <- env1[nbhd_index] # Make env1/nbhd indexing consistent
        names(env1) <- seq_along(nbhd_index)
        sim_steps <- 1
        objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
        par_true <- par0
        if (any(is.na(par_true))) return(NA)
        ll1 <- log_likelihood(par_true, objects1, debug = TRUE)
        ll2 <- log_likelihood(unlist(par_fitted[i,]), objects1, debug = TRUE)
        ll3 <- log_likelihood(c(0, 0, 0), objects1, debug = TRUE)

        l1 <- ll1[[2]][5,]
        l2 <- ll2[[2]][5,]
        l3 <- ll3[[2]][5,]

        hist(l2 - l1, 100, main = "Fitted - True", border = NA)
        abline(v = 0, lty = 2, lwd = 2)
    }
    )    
    out <- data.frame(id = seq_along(llike), ll_fit = llike, ll_0 = out)
    saveRDS(out, "optimtest.rds")


    # Why does leave one step out not work?
    step2 <- matrix(nrow = 5, ncol = 5, data = runif(25))
    step2 <- step2 / sum(step2)
    step1 <- step2[2:4, 2:4]
    step1 <- step1 / sum(step1)
    
    firsts <- sample(1:9, 800, replace = T, prob = as.vector(step1))
    emp1 <- table(firsts) / sum(table(firsts))
    plot(as.vector(step1), emp1)
    abline(0, 1)

    map <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)

    emp2 <- sapply(firsts, function(x) {
        start <- map[x]
        ind <- c(-6, -5, -4, -1, 0, 1, 4, 5, 6) + start
        val <- step2[ind] / sum(step2[ind])

    })

}

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
                                dir = simdir))
    ll2 <- -unlist(load_if_exists(paste0("ll_fit2", 1:sim_n, ".rds"), 
                                dir = simdir))
    aic1 <- 2 * 3 - 2 * ll1
    aic2 <- 2 * 3 - 2 * ll2
    # CHECK NUMBER OF PARAMETERS ^
    par(mfrow = c(1, 1))
    plot(aic1, aic2)
    abline(0, 1)
}