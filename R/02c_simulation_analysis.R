source("R/00_functions.R")

## Switches ====================================================================
plot_aic    <- F
param_plots <- T
simdir      <- "simulations/s4/"
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
    lines(x1, yhat2, lwd = 3, col = "red")
    for (i in seq_len(nrow(fit))) {
        lines(x1, y1[[i]], col = rgb(0, 0, 1, 0.1), lwd = 1)
        # lines(x1, y2[[i]], col = rgb(1, 0, 0, 0.5), lwd = 1.5)
        env <- points[[i]]$env
        move <- points[[i]]$move
        mw <- moving_window(env, move, window = 0.5)
        lines(mw$x, mw$y, col = rgb(0, 0, 0, 0.2), lwd = 1)
        abline(h = exp01(params$par_move), lty = 2, col = "red")
        model <- glm(move ~ env, family = binomial)
    }
}

mm <- unlist(lapply(paths, function(x) mean(x$att, na.rm = T)))


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