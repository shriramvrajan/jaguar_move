  
  ## Attraction function 1: environmental variables + home range ---------------
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6] +
  #                  par[7] * env$home)
  # attract_h <- exp(par[7] * env$home)
  # attract <- normalize_nbhd(attract_e * attract_h) #* normalize_nbhd(attract_t)
  # attract <- normalize_nbhd(attract_e)

  ## Attraction function 2: just home range ------------------------------------
  # attract_h <- exp(par[1] * env$home)
  # attract <- normalize_nbhd(attract_h)   
  
  # Attraction function 3: simulations with move param -------------------------
  # attract1 <- normalize_nbhd(exp(par[1] * env)) # + exp(par[2] * env2)
  # print(par)
  # attract1 <- env_function(env, par[2:4])
  # move_prob <- exp01(par[1])
  # attract <- t(apply(attract1, 1, function(r) {
  #   cent <- ceiling(length(r) / 2)
  #   r[cent] <- r[cent] * (1 - move_prob)
  #   r[-cent] <- r[-cent] * ((move_prob) / (sum(!is.na(r)) - 1))
  #   return(r / sum(r, na.rm = TRUE))
  # }))

  # Attraction function 4: With 0-1 parameter ----------------------------------
  # move_prob <- exp01(par[7])
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])
  # # attract_h <- exp(par[8] * env$home)
  # # attract <- normalize_nbhd(attract_e * attract_h)
  # attract <- normalize_nbhd(attract_e)
  # attract <- t(apply(attract, 1, function(r) {
  #   cent <- ceiling(length(r) / 2)
  #   r[cent] <- r[cent] * (1 - move_prob)
  #   r[-cent] <- r[-cent] * (move_prob / (sum(!is.na(r)) - 1))
  #   return(r / sum(r, na.rm = TRUE))
  # }))


make_movement_kernel <- function(n = 10000, sl_emp, ta_emp, max_dist, 
                                 minimum = 0, scale = 1000) {
  
  avail <- data.frame(sl = sample(sl_emp, n, replace = TRUE),
                      ta = sample(ta_emp, n, replace = TRUE))
  avail$x <- avail$sl * cos(avail$ta) / scale # scales from meters to km
  avail$y <- avail$sl * sin(avail$ta) / scale 
  avail$xi <- sapply(avail$x, function(x) ifelse(x > 0, floor(x), ceiling(x)))
  avail$yi <- sapply(avail$y, function(y) ifelse(y > 0, floor(y), ceiling(y)))
  density <- tapply(avail$x, list(avail$yi, avail$xi), length)
  density <- reshape2::melt(density)
  names(density) <- c("x", "y", "n")
  size <- max_dist * 2 + 1
  out <- data.frame(x = rep(-max_dist:max_dist, times = size),
                    y = rep(-max_dist:max_dist, each = size))

  out <- merge(out, density, by = c("x", "y"), all.x = TRUE)
  out$n[is.na(out$n)] <- 0
  out$n <- out$n + minimum
  out$n <- out$n / sum(out$n)
  
  return(out$n)
}

make_movement_kernel1 <- function(n = 10000, sl_emp, ta_emp, max_dist, minimum = 0) {
  
  # Fit gamma parameters? 
  avail <- data.frame(sl = sample(sl_emp, n, replace = TRUE),
                      ta = sample(ta_emp, n, replace = TRUE))
  avail$x <- avail$sl * cos(avail$ta) / 1000
  avail$y <- avail$sl * sin(avail$ta) / 1000
  avail$d <- floor(sqrt(avail$x^2 + avail$y^2))

  stay_prob <- length(which(avail$d == 0)) / length(avail$d)  

  moves <- avail$d[-which(avail$d == 0)]
  rate <- fitdist(moves, distr = "exp", method = "mle")$estimate
  # Return out$n as a vector of probabilities
  # But using 0-1 for the center cell + gamma dist for the rest

  # OK YOU CAN DO THIS GANBATTE !!!!!!!!!!!!!

  out <- data.frame(x = rep(-max_dist:max_dist, times = size),
                    y = rep(-max_dist:max_dist, each = size))
  return(out$n)

  # avail$xi <- sapply(avail$x, function(x) ifelse(x > 0, floor(x), ceiling(x)))
  # avail$yi <- sapply(avail$y, function(y) ifelse(y > 0, floor(y), ceiling(y)))
  # density <- tapply(avail$x, list(avail$yi, avail$xi), length)
  # density <- reshape2::melt(density)
  # names(density) <- c("x", "y", "n")
  # size <- max_dist * 2 + 1

  # out <- merge(out, density, by = c("x", "y"), all.x = TRUE)
  # out$n[is.na(out$n)] <- 0
  # out$n <- out$n + minimum
  # out$n <- out$n / sum(out$n)  
}

# Return -max(log(likelihood)) given parameters and model objects
# traditional SSF (log_likelihood0) and all others (log_likelihood)
log_likelihood0 <- function(par, objects) {
  # par        : Initial values of parameters for optimization
  env <- objects[[1]]
  max_dist <- objects[[2]]
  mk <- objects[[3]]
  obs <- objects[[4]]
  n_obs <- length(obs) + 1
  
  # Attractiveness function 0: traditional SSF 
  # attract_e <- exp(par[1] * env[, 1] + par[2] * env[, 2] + par[3] * env[, 3] +
  #                  par[4] * env[, 4] + par[5] * env[, 5] + par[6] * env[, 6])

  # Attractiveness function 1: simulations
  attract_e <- normalize_nbhd(env_function(env, par[2:4]))
  ncell_local <- (max_dist * 2 + 1) ^ 2
  
  p_obs <- sapply(seq_len(n_obs - 1), function(t) {
    env_local <- attract_e[(ncell_local * (t - 1) + 1):(ncell_local * t)]
    env_local <- env_local / sum(env_local, na.rm = TRUE)
    # TESTING WITHOUT MK 
    # p <- env_local
    # With MK 
    env_weight <- exp01(par[1])
    p <- env_weight * env_local + (1 - env_weight) * mk # mk = movement kernel
    return(ifelse((p[obs[t]] == 0), 0.0001, p[obs[t]]))
  })
  ll <- -sum(log(p_obs))
  if (is.infinite(ll)) ll <- 0
  # print(ll)
  return(ll)
}


if (refit_model0) {
  # Traditional SSF for comparison
  track <- make_full_track(id)
  sl_emp <- as.vector(na.exclude(track$sl))
  ta_emp <- as.vector(na.exclude(track$ta))
  mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, max_dist = max_dist)
}