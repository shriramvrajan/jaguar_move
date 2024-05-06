if (!exists("loaded")) source("R/00_functions.R") else print("R/00 already loaded")
loaded <- TRUE

## Switches ====================================================================
plot_aic       <- F
param_plots    <- F
debug_fit      <- F
indiv_analysis <- F

indiv          <- 7
simdir         <- "simulations/s1/"

parallel_setup(1)

## Load data ===================================================================
message("Loading data")
params       <- readRDS(paste0(simdir, "params.rds"))   
sim_interval <- params$sim_interval
sim_n        <- params$sim_n

paths <- readRDS(paste0(simdir, "paths.rds"))
jag_traject_cells <- lapply(paths, function(p) {
    tr <- cbind(p$x, p$y)
    ind <- seq(1, nrow(tr), sim_interval + 1)
    tr <- tr[ind, ]
    return(raster::cellFromXY(env01[[1]], tr[, 1:2]))
})

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

# dist <- lapply(jag_traject, function(tr) {
#     out <- c(0, sqrt(diff(tr[, 1])^2 + diff(tr[, 2])^2))
#     return(out)
# })
# max_dist <- ceiling(max(unlist(dist)) * 1.5)

## Individual analysis =========================================================

if (indiv_analysis) {

    i <- indiv
    max_dist <- 1
    step_range <- (2 * max_dist + 1) ^ 2
    buffersize <- 1

    ## Find probabilities for TRUE PAR and FITTED PAR and see what happens
    par_true <- unlist(params[10:12]) %>% as.numeric
    par_fit <- as.numeric(fit[i, 1:3])

    traject <- jag_traject_cells[[i]]
    prep_model_objects(traject, max_dist, env01[[1]])
 
    sim_steps <- sim_interval + 2
    objects1 <- list(env, nbhd, max_dist, sim_steps, to_dest, obs)

    l1 <- log_likelihood(par_true, objects1, debug = TRUE)
    l2 <- log_likelihood(par_fit, objects1, debug = TRUE)

    ### COMPARE LOG LIKELIHOOD LOGIC

    path_i <- paths[[i]]
    path_i$cell[1] <- (path_i$x[1] - 1) * 400 + path_i$y[1]
    nbhd_i <- make_nbhd(i = path_i$cell, sz = buffersize, r = env01[[1]], 
                        rdf = env01[[2]])
    env_i <- matrix(nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
    env_i[] <- env01[[2]]$sim1[nbhd_i]

    # nbhd <- nbhd_i
    att_true <- 1 / (1 + exp(par_true[1] + par_true[2] * env_i + 
                             par_true[3] * env_i^2))
    att_true <- att_true / rowSums(att_true)
    att_fit <- 1 / (1 + exp(par_fit[1] + par_fit[2] * env_i + 
                            par_fit[3] * env_i^2))
    att_fit <- att_fit / rowSums(att_fit)

    map <- 1:9
    names(map) <- c(-401, -400, -399, -1, 0, 1, 399, 400, 401)
    obs_ind <- path_i$cell[2:nrow(path_i)] - path_i$cell[1:(nrow(path_i) - 1)]
    obs_ind <- map[as.character(obs_ind)]

    p_true <- sapply(1:length(obs_ind), function(i) att_true[i, obs_ind[i]])
    p_fit <- sapply(1:length(obs_ind), function(i) att_fit[i, obs_ind[i]])

    ## Why don't these match above ones???  DEBUG LOGLIKELIHOOD 
    print(paste("Log-likelihood for true parameters: ", sum(log(p_true))))
    print(paste("Log-likelihood for fitted parameters: ", sum(log(p_fit))))

    debug(log_likelihood)

    par(mfrow = c(1, 1))
    plot(env_i, att_true, pch = 19, col = rgb(0, 0, 0, 1), 
         cex = 0.5, ylim = c(0, 0.25))
    points(env_i, att_fit, pch = 19, col = rgb(0, 0, 1, 0.2), cex = 0.5)
    abline(h = 1/9, col = "red")


}


## Explore parameter values ====================================================
if (param_plots) {
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

    points <- lapply(fit$id, function(i) {
        path <- paths[[i]]
        path$move <- c(0, sqrt(diff(path$x)^2 + diff(path$y)^2))
        path$move[path$move > 0] <- 1
        path$env <- env01[[2]]$sim1[cellFromXY(env01[[1]], path[, 1:2])]
        return(path[, c("move", "env")])
    })
    # Generated + fitted, all on same plot
    par(mfrow = c(1, 1))
    plot(x1, y0, type = "l", lwd = 3, ylim = c(0, 1))
    lines(x1, yhat, lwd = 3, col = "blue")
    lines(x1, yhat2, lwd = 3, col = rgb(1, 0, 0, 0.8))
    for (i in seq_len(nrow(fit))) {
        lines(x1, y1[[i]], col = rgb(0, 0, 1, 0.3), lwd = 1)
        readline("Press [enter] to continue")
    }

    par(mfrow = c(1, 3))
    plot(fit$c1)
    abline(h = par0[1], col = "red")
    plot(fit$b1)
    abline(h = par0[2], col = "red")
    plot(fit$a1)
    abline(h = par0[3], col = "red")
}

## Explore curve fitting =======================================================

logistic_plot <- function(par) {
    x1 <- seq(0, 8, length.out = 100)
    y0 <- 1 / (1 + exp(par0[1] + par0[2] * x1 + par0[3] * x1^2)) 
    y1 <- 1 / (1 + exp(par[1] + par[2] * x1 + par[3] * x1^2))
    plot(x1, y0, type = "l", lwd = 3, ylim = c(0, 1))
    lines(x1, y1, lwd = 3, col = "blue")
}

plot_samples <- function(n, par = par0) {
    x1 <- seq(0, 8, length.out = 10000)
    y1 <- 1 / (1 + exp(par[1] + par[2] * x1 + par[3] * x1^2))
    plot(x1, y1, col = "white", ylim = c(0, 1))
    # plot(x1, y1, type = "l", col = "blue", lwd = 3)
    for (i in seq_len(n)) {
        ind <- sort(sample(seq_len(length(x1)), 9))
        yy <- y1[ind] / sum(y1[ind])
        lines(x1[ind], yy, lwd = 1)
    }
}