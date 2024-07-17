  
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