
if (!exists("simname")) source("R/00_functions.R") else print("R/00 already loaded")
loaded <- TRUE

## Switches ====================================================================
plot_aic       <- F
param_plots    <- T
debug_fit      <- F
indiv_analysis <- F

simname        <- "s1"
simdir         <- paste0("simulations/", simname, "/")

## Load data ===================================================================
message("Loading data")
params       <- readRDS(paste0(simdir, "params.rds"))   
obs_interval <- params$obs_interval
sim_n        <- params$sim_n
step_size    <- params$step_size

paths <- readRDS(paste0(simdir, "paths.rds"))

env01 <- rast(paste0(simdir, "env01.tif"))
print(params)

fit1 <- load_if_exists(paste0("par_out1_", 1:sim_n, ".rds"), dir = simdir) %>%
        do.call(rbind, .) %>% 
        as.data.frame() 
fit2 <- load_if_exists(paste0("par_out2_", 1:sim_n, ".rds"), dir = simdir) %>%
        do.call(rbind, .) %>% 
        as.data.frame() 
fit <- cbind(fit1, fit2)
# if (all(is.na(fit))) fit <- readRDS(paste0(simdir, "par_out_all.rds"))
fit$id <- seq_len(nrow(fit))
names(fit) <- c("c1", "b1", "a1", "c2", "b2", "a2", "id")

if (any(is.na(fit$c1))) {
    posna <- which(is.na(fit$c1)) # positions of NA
    fit <- fit[-posna, ]
}

## Explore parameter values ====================================================
if (param_plots) {

    jag_traject <- lapply(paths, function(p) {
        out <- cbind(p$x, p$y)
        ind <- seq(1, nrow(out), obs_interval + 1)
        out <- out[ind, ]
        return(out)
    })
    jag_traject_cells <- lapply(jag_traject, function(tr) {
        out <- terra::cellFromRowCol(env01, tr[, 1], tr[, 2])
        return(out)
    })
    max_dist     <- step_size * (obs_interval + 1)

    message("Building global neighborhood")
    nbhd0 <- make_nbhd(i = seq_len(ncell(env01)), sz = step_size, r = env01,
                       rdf = raster_to_df(env01)) 

    objects_all <- lapply(jag_traject_cells, function(traject) {
        return(prep_model_objects(traject, max_dist, raster_to_df(env01), sim = TRUE))
    })

    par0 <- c(3, -2, 0.3)
    y0 <- plot_curve(par0, values = TRUE)    
    y1 <- lapply(seq_len(nrow(fit)), function(i) {
        mu <- objects_all[[i]]$mu_env
        sd <- objects_all[[i]]$sd_env
        out <- plot_curve(unlist(fit[i, 1:3]), mu = mu, sd = sd, values = TRUE)
        return(out)
    })

    y2 <- lapply(seq_len(nrow(fit)), function(i) {
        mu <- objects_all[[i]]$mu_env
        sd <- objects_all[[i]]$sd_env
        out <- plot_curve(unlist(fit[i, 4:6]), mu = mu, sd = sd, values = TRUE)
        return(out)
    })

    # Generated + fitted, all on same plot
    plotpdf(nm = paste0("figs/sims/", simname, "plot.pdf"), x = 8, y = 4)
    par(mfrow = c(1, 2))
    plot(y0, type = "l", lwd = 3, ylim = c(0, 1), xlab = "Environmental variable",
         ylab = "Attraction", main = "Step selection")
    # lines(x1, yhat, lwd = 3, col = "#1a1a9e")
    # lines(x1, yhat2, lwd = 3, col = rgb(1, 0, 0, 0.8))
    for (i in seq_len(nrow(fit))) {
        lines(y1[[i]], col = rgb(0, 0, 1, 0.3), lwd = 3)
        # readline(paste(i, "Press [enter] to continue"))
    }
    
    plot(y0, type = "l", lwd = 3, ylim = c(0, 1), xlab = "Environmental variable",
         ylab = "Attraction", main = "Path propagation")
    for (i in seq_len(nrow(fit))) {
        lines(y2[[i]], col = rgb(0, 0, 1, 0.3), lwd = 3)
        # readline(paste(i, "Press [enter] to continue"))
    }
    dev.off()

    plotpdf(nm = "figs/simplot2.pdf", x = 12, y = 4)
    par(mfrow = c(1, 3))
    plot(fit$c1)
    abline(h = par0[1], col = "red")
    plot(fit$b1)
    abline(h = par0[2], col = "red")
    plot(fit$a1)
    abline(h = par0[3], col = "red")
    dev.off()
}

## Individual analysis =========================================================

if (indiv_analysis) {
    par(mfrow = c(1, 1))
    
    jag_traject_cells <- lapply(paths, function(p) {
        tr <- cbind(p$x, p$y)
        ind <- seq(1, nrow(tr), obs_interval + 1)
        tr <- tr[ind, ]
        return(raster::cellFromRowCol(env01, tr[, 1], tr[, 2]))
    })    

    max_dist     <- step_size * (obs_interval + 1) # max dist between obs
    sim_steps    <- max_dist + 1                   # number steps to simulate
    ncell_local  <- (2 * max_dist + 1) ^ 2         # local nbhd size

    comp_hist <- function(i) {
        par_true <- unlist(params[10:12]) %>% as.numeric # true parameters
        par_fit <- as.numeric(fit[i, 1:3])               # fitted parameters
        traject <- jag_traject_cells[[i]]                # individual trajectory
        
        prep_model_objects(traject, max_dist, env01) # prepare model objects
        objects1 <- list(env, nbhd, max_dist, sim_steps, to_dest, obs)

        l1 <- log_likelihood(par_true, objects1, debug = TRUE)
        l2 <- log_likelihood(par_fit, objects1, debug = TRUE)
        print(paste("Log-likelihood for true parameters: ", l1[[1]]))
        print(paste("Log-likelihood for fitted parameters: ", l2[[1]]))

        # Plot env vs probability
        par(mfrow = c(1, 2))
        probs <- paths[[i]]$att[-1]
        names(probs) <- paths[[i]]$cell[-1]
        x <- env01[traject][[1]]
        y0 <- probs[traject]
        y1 <- l1[[2]][sim_steps, ]
        y2 <- l2[[2]][sim_steps, ]
        plot(x, y1)
        abline(h = 1 / 9, lty = 2, lwd = 3, col = "red")
        plot(x, y2)
        abline(h = 1 / 9, lty = 2, lwd = 3, col = "red")

        # Plot histograms
        range1 <- range(c(l1[[2]][sim_steps, ], l2[[2]][sim_steps, ]), na.rm = T)
        hist(l2[[2]][sim_steps, ], col = rgb(1, 0, 0, 0.5), 50, border = NA, 
            xlim = range1, main = paste(i, "Red = fitted, blue = true"))
        hist(l1[[2]][sim_steps, ], col = rgb(0, 0, 1, 0.5), 100, border = NA, add = T)

        chance <- 1 / ((sim_steps - 1) * 2 + 1) ^ 2
        abline(v = chance, lty = 2, lwd = 3)
    }

    comp_hist(2)
}

## Parameter landscape =========================================================

if (FALSE) {
jag_traject <- lapply(paths, function(p) {
    out <- cbind(p$x, p$y)
    ind <- seq(1, nrow(out), obs_interval + 1)
    out <- out[ind, ]
    return(out)
})
jag_traject_cells <- lapply(jag_traject, function(tr) {
    out <- terra::cellFromRowCol(env01, tr[, 1], tr[, 2])
    return(out)
})
# Make global neighborhood from raster
message("Building global neighborhood")
nbhd0 <- make_nbhd(i = seq_len(ncell(env01)), sz = step_size, r = env01) 

max_dist <- step_size * (obs_interval + 1)
ncell_local <- (2 * max_dist + 1) ^ 2
sim_steps   <- obs_interval * step_size + 2
# Number of steps to simulate, interval + first and last steps

objects_all <- lapply(jag_traject_cells, function(traject) {
    return(prep_model_objects(traject, max_dist, env01))
})

ll_product <- function(par) {
   ll <- lapply(objects_all, function(o) log_likelihood(par, objects = o)) %>%
            unlist() %>%
            sum()
}

par_test <- expand.grid(c1 = 2,
                        b1 = seq(-5, 5, length.out = 20),
                        a1 = seq(-5, 5, length.out = 20))

ll_test <- apply(par_test, 1, function(p) {
    print(p)
    ll_product(p) %>% return()
})

p <- ggplot(data.frame(ll = ll_test, par_test), aes(x = b1, y = a1, fill = ll)) +
    geom_tile() +
    scale_fill_viridis_c() +
    theme_minimal() +
    labs(title = "Log-likelihood landscape", x = "b1", y = "a1")
}



### Scratch ===================================================================

## Reconstructing the log-likelihood function
# path_i <- paths[[i]]
# path_i$cell[1] <- (path_i$x[1] - 1) * 400 + path_i$y[1]
# nbhd_i <- make_nbhd(i = path_i$cell, sz = step_size, r = env01)
# env_i <- matrix(nrow = nrow(nbhd_i), ncol = ncol(nbhd_i))
# env_i[] <- env01[nbhd_i]

# # nbhd <- nbhd_i
# att_true <- 1 / (1 + exp(par_true[1] + par_true[2] * env_i + 
#                          par_true[3] * env_i^2))
# att_true <- att_true / rowSums(att_true)
# att_fit <- 1 / (1 + exp(par_fit[1] + par_fit[2] * env_i + 
#                         par_fit[3] * env_i^2))
# att_fit <- att_fit / rowSums(att_fit)

# map <- 1:9
# names(map) <- c(-401, -400, -399, -1, 0, 1, 399, 400, 401)
# obs_ind <- path_i$cell[2:nrow(path_i)] - path_i$cell[1:(nrow(path_i) - 1)]
# obs_ind <- map[as.character(obs_ind)]

# p_true <- sapply(1:length(obs_ind), function(i) att_true[i, obs_ind[i]])
# p_fit <- sapply(1:length(obs_ind), function(i) att_fit[i, obs_ind[i]])

# ## Why don't these match above ones???  DEBUG LOGLIKELIHOOD 
# print(paste("Log-likelihood for true parameters: ", sum(log(p_true))))
# print(paste("Log-likelihood for fitted parameters: ", sum(log(p_fit))))

# par(mfrow = c(1, 1))
# plot(env_i, att_true, pch = 19, col = rgb(0, 0, 0, 1), 
#      cex = 0.5, ylim = c(0, 0.25))
# points(env_i, att_fit, pch = 19, col = rgb(0, 0, 1, 0.2), cex = 0.5)
# abline(h = 1 / 9, col = "red")

### Explore curve fitting

# logistic_plot <- function(par) {
#     x1 <- seq(0, 8, length.out = 100)
#     y0 <- 1 / (1 + exp(par0[1] + par0[2] * x1 + par0[3] * x1^2)) 
#     y1 <- 1 / (1 + exp(par[1] + par[2] * x1 + par[3] * x1^2))
#     plot(x1, y0, type = "l", lwd = 3, ylim = c(0, 1))
#     lines(x1, y1, lwd = 3, col = "blue")
# }

# plot_samples <- function(n, par = par0) {
#     x1 <- seq(0, 8, length.out = 10000)
#     y1 <- 1 / (1 + exp(par[1] + par[2] * x1 + par[3] * x1^2))
#     plot(x1, y1, col = "white", ylim = c(0, 1))
#     # plot(x1, y1, type = "l", col = "blue", lwd = 3)
#     for (i in seq_len(n)) {
#         ind <- sort(sample(seq_len(length(x1)), 9))
#         yy <- y1[ind] / sum(y1[ind])
#         lines(x1[ind], yy, lwd = 1)
#     }
# }
