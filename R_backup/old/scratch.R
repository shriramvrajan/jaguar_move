
#   # ## Extra variables
#   # df$nmove <- sapply(1:njag, function(x) {
#   #   length(which(jag_move$ID == as.numeric(jag_id[x])))
#   # })
#   # df$ndays <- sapply(1:njag, function(x) {
#   #   moves <- jag_move[ID == as.numeric(jag_id[x])]
#   #   dates <- sort(as.Date(sapply(moves$timestamp, function(dt) {
#   #           strsplit(as.character(dt), " ")[[1]][1]
#   #       }), format = "%m/%d/%y"))
#   #   return(as.numeric(difftime(dates[length(dates)], dates[1])))
#   # })
#   # df$meandist <- tapply(jag_move$dist, jag_move$ID, function(x) mean(x, na.rm = TRUE))
#   # df$totdist <- tapply(jag_move$dist, jag_move$ID, function(x) sum(x, na.rm = TRUE))

# # pm <- 0.7
# # step_size <- 1

# #   aa1 <- t(apply(att, 1, function(env) {
# #     cent    <- ceiling(length(env) / 2)
# #     missing <- length(is.na(env))
# #     if (missing > 0) {
# #       kernel[-cent] <- (1 - kernel[cent]) / (length(kernel) - missing - 1)
# #     }
# #     p <- env * kernel
# #     return(p / sum(p, na.rm = T))    
# #   }))

# # # p = 1 / (1 + 8 * exp(-alpha * 1))
# # alpha <- -log((1/pm - 1)/8)


# #   kernel0 <- dexp(0:step_size, rate = alpha) # move par
# #   kernel <- matrix(0, nrow = step_size * 2 + 1, ncol = step_size * 2 + 1)
# #   center <- step_size + 1
# #   for (i in seq_len(center + step_size)) {
# #     for (j in seq_len(center + step_size)) {
# #       kernel[i, j] <- kernel0[max(c(abs(i - center), abs(j - center))) + 1]
# #     }
# #   }
# #   kernel <- kernel / sum(kernel)

# #     # cent <- ceiling(length(env) / 2)
# #     # env[cent] <- env[cent] * (1 - move_prob)
# #     # env[-cent] <- env[-cent] * (move_prob / (sum(!is.na(env)) - 1))
# #     # return(env / sum(env))

# # Plotting from results ------------------------------------------------------
# # if (plot_res) {
# #     # vars: 1 footprint, 2 elev, 3 slope, 4 forestcover, 5 distwater, 
# #     #       6 distroad, 7 homerange
# #     evars <- c("footprint", "elev", "slope", "forest", "distance to water", "distance to road",
# #             "unfamiliarity")

# #     i <- 3
    
# #     barplot(tapply(res[, paste0("p", i)], res$bio, function(x) mean(x, na.rm = TRUE)),
# #             col = "#a37373", border = NA, main = evars[i])



# #     res <- as.data.frame(res) # still not used to plotting with tibbles

    

# #     par(mfrow = c(2, 4))

# #     for (i in seq_len(length(evars))) {
# #         hist(res[, paste0("p", i)], 20, col = "#a37373", border = NA,
# #              xlab = "Parameter value", main = evars[i])
# #         abline(v = 0, lwd = 2, col = "black")
# #         print(evars[i])
# #         print(summary(res[, paste0("p", i)]))
# #     }
# # }

# ## Markov version of jag_path --------------------------------------------------

# # Generate a jaguar path of n steps starting from (x0, y0) with environmental
# # preference parameters par[] and search neighborhood size neighb
# jag_path <- function(x0, y0, n_step, par, neighb, type = 1, tprob = 0) {
#     # type: 1 = env1 only, 2 = multi-state model

#     # Set initial state
#     path <- matrix(NA, nrow = n_step, ncol = switch(type, 4, 6))
#     path[1, ] <- switch(type, 
#                         c(x0, y0, NA, NA), c(x0, y0, NA, state0, NA, NA))

#     # Probability of moving from current grid cell
#     move_prob <- exp(par[2]) / (1 + exp(par[2]))
    
#     if (type == 2) {
#         # Transition probabilities: p12, p21, p11, p22
#         tprob <- c(tprob, 1 - tprob)
#     }

#     for (i in 2:n_step) {
#         if (i %% 10 == 0) print(i)
#         pos <- path[i - 1, 1:2]
#         if (type == 2 && i %% sim_interval == 0) {
#             state <- path[i - 1, 4]
#             # Transition to new state
#             if (state == 1) {
#                 if (runif(1) < tprob[1]) state <- 2
#             } else {
#                 if (runif(1) < tprob[2]) state <- 1
#             }
#         }
#         nbhd <- as.vector(make_nbhd(r = env01[[1]], rdf = env01[[2]], sz = neighb,
#                           i = cellFromRowCol(env01[[1]], pos[1], pos[2])))
#         a1 <- sapply(nbhd, function(x) {
#           if (is.na(x)) return(NA) else return(exp(env01[[1]][x] * par[1])$sim1)
#         })
#         if (any(is.na(a1))) a1[is.na(a1)] <- 0
#         a1 <- a1 / sum(a1)
#         cent <- ceiling(length(a1) / 2)
#         a1[cent] <- a1[cent] * (1 - move_prob)
#         a1[-cent] <- a1[-cent] * (move_prob / (sum(!is.na(a1)) - 1))
#         attract <- a1 / sum(a1)
        
#         if (type == 2) {
#           a2 <- exp(env02[[1]][nbhd] * par[2])$sim1
#           if (any(is.na(a2))) a2[is.na(a2)] <- 0
#           a2 <- a2 / sum(a2)

#           attract <- switch(state, a1, a2)
#         }

#         step <- sample(seq_len(length(attract)), 1, prob = attract)
#         path[i, ] <- switch(type,
#                             c(rowColFromCell(env01[[1]], nbhd[step]), 
#                               nbhd[step], a1[step]),
#                             c(rowColFromCell(env01[[1]], nbhd[step]), 
#                               nbhd[step], state, a1[step], a2[step]))
#     }
#     path <- as.data.frame(path)
#     names(path) <- switch(type, c("x", "y", "cell", "a1"), 
#                           c("x", "y", "cell", "state", "a1", "a2"))
#     return(path)
# }

# # Two state simulation =========================================================

# probs <- load_if_exists(paste0("p", 1:sim_n, ".RDS"), 
#                         dir = "data/output/simulations")
# currents   <- load_if_exists(paste0("current", 1:sim_n, ".RDS"), 
#                         dir = "data/output/simulations")
# par(mfrow = c(1, 3))
# res <- lapply(1:sim_n, function(i) {
#     path <- paths[[i]]
#     prob <- probs[[i]]
#     curr1  <- currents[[i]][[1]]
#     curr2  <- currents[[i]][[2]]
#     side <- sqrt(dim(curr1)[1])

#     state <- path$state[seq(1, nrow(path), sim_interval)]
#     p     <- prob[nrow(prob), ]

#     s <- 99
#     ss <- 10
#     r1 <- raster(nrows = side, ncols = side, vals = curr1[, s, ss])
#     r2 <- raster(nrows = side, ncols = side, vals = curr2[, s, ss])
#     terra::plot(r1, main = paste0("Kernel #", i))
#     terra::plot(r2, main = paste0("Kernel #", i))
#     hist(p[state == 1], 11, main = paste0("Path #", i), col = rgb(0, 0, 1, 0.5),
#          xlim = c(0, max(p, na.rm = TRUE)), border = NA)
#     hist(p[state == 2], 11, add = TRUE, col = rgb(1, 0, 0, 0.5), border = NA)
#     abline(v = mean(p[state == 1], na.rm = TRUE), col = "blue")
#     abline(v = mean(p[state == 2], na.rm = TRUE), col = "red")
#     abline(v = 1 / (21 ^ 2))

#     # plot_path(path)

#     # x <- readline("Press [enter] to continue") 

#     # a1_t <- path$a1[which(path$state == 1)]
#     # a1_f <- path$a2[which(path$state == 1)]
#     # a2_t <- path$a2[which(path$state == 2)]
#     # a2_f <- path$a1[which(path$state == 2)]

#     # breaks1 <- (0:25 * 0.04) * max(c(a1_t, a1_f), na.rm = TRUE)
#     # hist(a1_t, 30, col = rgb(0, 0, 1, 0.4), border = NA, breaks = breaks1)
#     # abline(v = mean(a1_t, na.rm = TRUE), col = "blue")
#     # hist(a1_f, 30, col = rgb(1, 0, 0, 0.4), add = TRUE, border = NA, 
#     # breaks = breaks1)
#     # abline(v = mean(a1_f, na.rm = TRUE), col = "red")

#     # breaks2 <- (0:25 * 0.04) * max(c(a2_t, a2_f), na.rm = TRUE)
#     # hist(a2_f, 30, col = rgb(1, 0, 0, 0.4), border = NA, breaks = breaks2)
#     # abline(v = mean(a2_f, na.rm = TRUE), col = "red")
#     # hist(a2_t, 30, col = rgb(0, 0, 1, 0.4), add = TRUE, border = NA,  
#     #      breaks = breaks2)
#     # abline(v = mean(a2_t, na.rm = TRUE), col = "blue")
# })

# ## Plotting of parameter landscapes (2D by color, only 2 pars) ==================

# gen_parscape  <- F
# plot_parscape <- F

# # Generate optim landscape 
# if (gen_parscape) {
#     message("Generating fitting landscape...")

#     x <- seq(-4, 4, length.out = 50)
#     y <- seq(-4, 4, length.out = 50)
#     testvals <- expand.grid(x, y)
#     names(testvals) <- c("x", "y")

#     to_gen <- seq_len(nrow(fit))[-which(is.na(fit$a1))]
#     for (i in to_gen) {
#     # foreach(i = to_gen) %dopar% {    
#         message(paste0("Generating landscape", i, " / ", nrow(fit)))
#         traject <- jag_traject_cells[[i]]
#         prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
#                         rdf = env01[[2]])
#         env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
#         # Normalizing desired environmental variables for extended neighborhood
#         env1 <- env1[nbhd_index]
#         # Make indexing consistent with env
#         names(env1) <- seq_len(length(nbhd_index))
#         sim_steps <- sim_interval * 2

#         objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)

#         message("Calculating log likelihoods")
#         testvals$ll <- unlist(foreach(j = seq_len(nrow(testvals)), 
#                               .export = c("dest")) %dopar% {
#             message(paste0("Fitting row #: ", j, " / ", nrow(testvals)))
#             val <- as.numeric(testvals[j, ])
#             ll <- log_likelihood(val, objects1)
#             message(paste0("Fitted row #: ", j, " / ", nrow(testvals)))
#             return(ll)
#         })
#         message("Aggregating results")
#         saveRDS(testvals, "data/output/simulations/testll.rds")
#     }
# }

# # Plot parameter landscape 
# if (plot_parscape) {
#         testvals <- readRDS("data/output/simulations/testll.rds")
#         plot <- ggplot(testvals, aes(x = x, y = y, fill = ll)) +
#             geom_tile() +
#             scale_fill_viridis_c(option = "turbo") +
#             theme_minimal() +
#             labs(title = paste0("Log likelihood surface", i),
#                 x = "a1", y = "b1") +
#             theme(legend.position = "none") +
#             geom_point(aes(x = 3, y = -2), color = "magenta", size = 3)
#         ggsave(paste0("data/output/simulations/ll_surface_", i, ".png"), plot,
#                bg = "white")
# }


# ## Debug =======================================================================

# if (debug_fit) {
#     # What was I trying to do here...
#     parallel_setup(ncore_fit)
#     llike <- load_if_exists(paste0(simdir, "ll_fit1", 1:sim_n, ".rds"), 
#                             dir = ".") %>%
#              unlist(.)
#     par_fitted <- load_if_exists(paste0(simdir, "par_out_", 1:sim_n, ".rds"), 
#                                  dir = ".") %>%
#                   do.call(rbind, .) %>%
#                   as.data.frame()
#     # out <- foreach(i = seq_along(llike), .combine = c) %dopar% {
#     out <- sapply(seq_along(llike), function(i) {
#         ## !! JUST TESTING WITH ONE FIRST, generalize later !!
#         i <- 12
#         traject <- jag_traject_cells[[i]]
#         prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
#                 rdf = env01[[2]])
#         env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))
#         # Normalizing desired environmental variables for extended neighborhood
#         env1 <- env1[nbhd_index] # Make env1/nbhd indexing consistent
#         names(env1) <- seq_along(nbhd_index)
#         sim_steps <- 1
#         objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
#         par_true <- par0
#         if (any(is.na(par_true))) return(NA)
#         ll1 <- log_likelihood(par_true, objects1, debug = TRUE)
#         ll2 <- log_likelihood(unlist(par_fitted[i,]), objects1, debug = TRUE)
#         ll3 <- log_likelihood(c(0, 0, 0), objects1, debug = TRUE)

#         l1 <- ll1[[2]][5,]
#         l2 <- ll2[[2]][5,]
#         l3 <- ll3[[2]][5,]

#         hist(l2 - l1, 100, main = "Fitted - True", border = NA)
#         abline(v = 0, lty = 2, lwd = 2)
#     }
#     )    
#     out <- data.frame(id = seq_along(llike), ll_fit = llike, ll_0 = out)
#     saveRDS(out, "optimtest.rds")
# }



# ## Plot model AIC ==============================================================
# if (plot_aic) {
#     par(mfrow = c(1, 2))
#     plot(fit$a1, fit$b1)
#     abline(h = 0)
#     abline(v = 0)
#     plot(fit$a2, fit$b2)
#     abline(h = 0)
#     abline(v = 0)

#     sim_n <- nrow(fit)
#     ll1 <- -unlist(load_if_exists(paste0("ll_fit1", 1:sim_n, ".rds"), 
#                                 dir = simdir))
#     ll2 <- -unlist(load_if_exists(paste0("ll_fit2", 1:sim_n, ".rds"), 
#                                 dir = simdir))
#     aic1 <- 2 * 3 - 2 * ll1
#     aic2 <- 2 * 3 - 2 * ll2
#     # CHECK NUMBER OF PARAMETERS ^
#     par(mfrow = c(1, 1))
#     plot(aic1, aic2)
#     abline(0, 1)
# }

# if (FALSE) {
#     ## Toy model for logic -----------------------------------------------------
#     # Global neighborhood:
#     step2 <- matrix(nrow = 5, ncol = 5, data = runif(25))
#     step2 <- step2 / sum(step2)
#     # Central cell is the starting point
#     step1 <- step2[2:4, 2:4]
#     # Local neighborhood for step 1:
#     step1 <- step1 / sum(step1)
#     # Sample for first step
#     firsts <- sample(1:9, 1000, replace = T, prob = as.vector(step1))
#     emp1 <- table(firsts) / sum(table(firsts))
#     plot(as.vector(step1), as.vector(emp1), main = "First step", 
#          xlab = "True", ylab = "Empirical", pch = 19)
#     abline(0, 1)
#     # Cell indices in global neighborhood for the central 3x3
#     map <- c(7, 8, 9, 12, 13, 14, 17, 18, 19)
#     # Sample for second step, given first step
#     emp2 <- sapply(firsts, function(x) {
#         start <- map[x]
#         ind <- c(-6, -5, -4, -1, 0, 1, 4, 5, 6) + start
#         val <- step2[ind] / sum(step2[ind])
#         step <- sample(1:9, 1, prob = val)
#         return(ind[step])
#     })
#     emp2 <- as.vector(table(factor(emp2, levels = 1:25)))
#     emp2 <- emp2 / sum(emp2)

#     # True probabilities for second step
#     joint <- lapply(map, function(i) {
#         ind <- c(-6, -5, -4, -1, 0, 1, 4, 5, 6) + i
#         val <- step2[ind] / sum(step2[ind])
#         val <- val * step1[which(map == i)]
#         names(val) <- ind
#         return(val)
#     })
#     true2 <- matrix(ncol = 25, nrow = 9)
#     for (i in seq_along(joint)) {
#         true2[i, as.numeric(names(joint[[i]]))] <- joint[[i]]
#     }
#     true2 <- colSums(true2, na.rm = T)
#     plot(as.vector(true2), emp2, main = "Second step", pch = 19,
#          xlab = "True", ylab = "Empirical")
#     abline(0, 1)

#     ## Joint probabilities?? ---------------------------------------------------
#     i <- 9
#     traject <- jag_traject_cells[[i]]
#     prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
#             rdf = env01[[2]])
#     env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))  
#     # Normalizing desired environmental variables for extended neighborhood
#     env1 <- env1[nbhd_index] # Make env1/nbhd indexing consistent
#     names(env1) <- seq_len(length(nbhd_index))
#     sim_steps <- sim_interval + 2
#     objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
    

#     opt0 <- optim(par = par0,       fn = log_likelihood, objects = objects1)
#     opt1 <- optim(par = c(1, 1, 1), fn = log_likelihood, objects = objects1)
    
#     print(opt0$par)
#     print(opt1$par)
    
#     ll0 <- log_likelihood(par0, objects1, debug1=T)
#     ll1 <- log_likelihood(opt0$par, objects1, debug1=T)
#     ll2 <- log_likelihood(opt1$par, objects1, debug1=T)

#     p0 <- ll0[[2]][2, ]
#     p1 <- ll1[[2]][2, ]

#     plot(fit$c1, main = "Fitted parameter values", pch = 19)
#     abline(h = 2, lty = 2, col = "red")

#     plot(p0, p1, xlab = "Step probabilities given fitted values", 
#          ylab = "Given true values", pch = 19)

#     refits <- lapply(seq_along(fit[, 1]), function(i) {
#         print(i)
#         current_jag <- i
#         traject <- jag_traject_cells[[i]]
#         prep_model_objects(traject, max_dist, nbhd0 = nbhd0, r = env01[[1]], 
#                 rdf = env01[[2]])
#         env1 <- scales::rescale(env01[[2]]$sim1[nbhd_index], to = c(0, 1))  
#         # Normalizing desired environmental variables for extended neighborhood
#         env1 <- env1[nbhd_index] # Make env1/nbhd indexing consistent
#         names(env1) <- seq_len(length(nbhd_index))
#         sim_steps <- sim_interval + 2
#         objects1 <- list(env1, nbhd, max_dist, sim_steps, to_dest, obs)
        
#         opt <- optim(par = c(1, 1, 1), fn = log_likelihood, objects = objects1)
#         return(opt$par)
#     }) %>% do.call(rbind, .) %>% as.data.frame()
#     names(refits) <- c("c1", "b1", "a1")

#     plot(fit$c1, refits$c1)
# }

# ### Stay/move logic

# Pm <- 0.8 / 8
# Ps <- 0.2

# env <- runif(9)
# env <- env / sum(env)
# moves <- c(1:4, 6:9)

# # first way to do it
# att1        <- env
# att1[moves] <- att1[moves] / sum(att1[moves]) * (1 - Ps)
# att1[5]     <- Ps

# # second way to do it - the current way
# att2 <- env * c(Pm, Pm, Pm, Pm, Ps, Pm, Pm, Pm, Pm)
# att2 <- att2 / sum(att2)


#   # move_prob <- exp01(par[length(par)]) # last one is move par
#   # attract <- t(apply(attract, 1, function(r) {
#   #   cent <- ceiling(length(r) / 2)
#   #   # r[cent] <- r[cent] * (1 - move_prob)
#   #   # r[-cent] <- r[-cent] * (move_prob / (sum(!is.na(r)) - 1))
#   #   # return(r / sum(r, na.rm = TRUE))
#   #   r[cent] <- move_prob
#   #   r[-cent] <- r[-cent] / sum(r[-cent], na.rm = TRUE) * (1 - move_prob)
#   #   return(r)
#   # })) # multiply environmental kernel and movement kernel


# ### Figures for committee report 2024 ==========================================

# plotpdf("figs/example.pdf", x = 8, y = 8)
# par(mfrow = c(2, 2))
# terra::plot(env01)
# plot_path(a)
# plot_curve(par0)
# dev.off()

# ## Scratch =====================================================================

# ## Traditional SSF fitting code
#     # message("Fitting parameters for model 2: traditional SSF") #------------ 
#     # obs_points <- as.data.frame(jag_traject[[i]]) 
#     # names(obs_points) <- c("x", "y")
#     # tr <- amt::steps(amt::make_track(obs_points, x, y))
#     # sl_emp <- as.vector(na.exclude(tr$sl_))
#     # ta_emp <- as.vector(na.exclude(tr$ta_))
#     # mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, 
#     #                            max_dist = max_dist, scale = 1)
#     # objects2 <- list(env, max_dist, mk, obs)
#     # par_out2 <- optim(par_start, log_likelihood0, objects = objects2)
#     # ll <- log_likelihood0(par_out2$par, objects2)
#     # message(paste0("Saving log-likelihood for model 2: ", i))
#     # saveRDS(ll, file = paste0("simulations/ll_fit2", i, ".rds"))