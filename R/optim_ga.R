# Genetic Algorithm integrated with optim to reduce the local minima problem by Brian Leung - March 2025.

library(doParallel)
library(foreach)

registerDoParallel()

wrapper <- function(par, le_func, df, control) {
    return(optim(par, le_func, df = df, control = control))
}

change <- function(par, p_change) {
    tp <- sample(1:3, ncore, replace = T, prob = p_change)
    for (i in 2:ncore) {
        if (tp[i] == 1) { # mutation
            par[i, ] <- par[i, ] * (.8 + .4 * runif(ncol(par)))
        } else if (tp[i] == 2) { # cross over
            s <- sample(1:ncore, 1) # find other parent
            t <- sample(seq_len(ncol(par)), 1) # choose 1 trait to swap
            par[i, t] <- par[s, t]
            t <- sample(seq_len(ncol(par)), 1) 
            # also mutate one other trait, so that make sure don't have exact duplicates
            par[i, t] <- par[i, t] * (.8 + .4 * runif(1))
        } else { # co-dominance
            s <- sample(1:ncore, 1) # find other parent
            par[i, ] <- (par[i, ] + par[s, ]) / 2
            t <- sample(seq_len(ncol(par)), 1) 
            # also mutate one other trait, so that make sure don't have exact duplicates
            par[i, t] <- par[i, t] * (.8 + .4 * runif(1))
        }
    }
    return(par)
}

optim_gm <- function(par, le_func, df = NULL, wrap_optim = wrapper, 
                     control = list(), ngen = 5, maxit = 50, ncore = 15, 
                     p_change = c(.5, .25, .25)) {
    # p_change - mutation,cross over, co-dominance 
    # 	operations of gm - create population, select who survives based on 
    # "fitness". Has mutation and cross-over (which in this case includes 
    # independent assortment). Can also mix the values of both (e.g, take midpoint)
    # keep the best performer. Choose the other ones based on their likelihoods.
    control$maxit <- maxit
    options(cores = ncore)

    for (gen in 1:ngen) {
        mem_optim <- foreach(i = 1:ncore) %dopar% {
                tmp_optim <- wrap_optim(par[i, ], le_func, df, control)
                return(unlist(c(lk = tmp_optim$val, tmp_optim$par)))
            }
        mem_optim <- as.data.frame(do.call(rbind, mem_optim))
        # now compare - keep the best one, allow the others to reproduce
        # 		par2=par
        min <- mem_optim[which.min(mem_optim$lk), ]
        par[1, ] <- min[, -1]
        # for each individual, decide which change or whether to do it
        # which individuals survive?
        p <- exp(-(mem_optim[, 1] - min$lk))
        s <- sample(1:ncore, ncore - 1, replace = T, prob = p) # these are the new individuals
        par[-1, ] <- mem_optim[s, -1] # take the parameter outcomes of each of the optims (remove the likelihood values)
        par <- change(par, p_change)
        print(Sys.time())
        print(c(gen, min$lk))
    }
    return(min)
}

# *** running the optim_gm
# **	set initialization parameters
ini_par <- data.frame(V1, V2)
# ** iterate through remaining cores allocated to project
# ** set DIFFERENT starting values
for (j in 2:ncore) {
    ini_par[j, ] <- c(ini_par[1, "V1"] * (1.5 * runif(1) + .5), 
                      ini_par[1, "V2"] * (1.5 * runif(1) + .5)) 
# the null model does not have b (which is the parameter of interest in our case)
}
# ** Note that if are comparing this against a previous model (e.g., simpler model), 
# ** can set the best optim run parameters as one of the starting values
# ** pass it parameter matrix to fit, the evaluation function, and a number of optional variables
# ** optional variables include the data, the wrapper function around optim, 
# ** control variables used in optim, number of generations, number of iterations 
# ** for optim, number of cores to use, and the probability of the different "evolutionary" operations
ret <- optim_gm(ini_par, le_func, df, wrapper, control, ngen, maxit)
