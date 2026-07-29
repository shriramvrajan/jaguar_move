source("R/functions.R")
source("R/classes.R")

test_warmstart <- FALSE

ssfn <- "data/output/emp_ss_1o_mm2.rds"
ppfn <- "data/output/emp_pp_1o_mm2.rds"
r1 <- results_set$new(r_ss = ssfn, r_pp = ppfn, env_type = "1o")$res_table

## test the forest-focal brdf

jag <- jaguar$new(63, results = r1)
sspar <- as.numeric(jag$results[paste0("ss_par", 1:9)])

brdf <- readRDS("data/env_layers.rds")
ssfit0 <- jag$fit_ss(par = sspar)

brdf <- readRDS("data/env_layers_fc9.rds")
ssfit1 <- jag$fit_ss(par = sspar)


### test whether warm starting with ss best fit is helping--------------------

ss_res <- readRDS(ssfn); pp_res <- readRDS(ppfn)
par0 <- ss_res[[which(jag_id$jag_id == id)]]$par
par1 <- pp_res[[which(jag_id$jag_id == id)]]$par

if (test_warmstart) {
    result    <- jag$fit_pp(par0, env_type = "1o")
    env_index <- head(seq_along(par0), -2) # -last 2
    par1 <- c(rnorm(length(env_index), 0, 0.1), tail(par0, 2))
    result2   <- jag$fit_pp(par1, env_type = "1o")
    saveRDS(list("warm" = result, "cold" = result2), "data/output/test.rds")
}

ssfit <- jag$fit_ss(par0)
ppfit <- jag$fit_pp(par0)

r1 <- jag$eval_ss(par = ssfit$par, env_type = "1o")
message(paste("SS log(L):", r1$ll_total))
r2 <- jag$eval_pp(par = ppfit$par, n_jump = ppfit$n_jump)
message(paste("PP log(L):", r2$ll_total, "| n_jump =", ppfit$n_jump,
              "| best_r =", r2$best_r))

## nesting check: multipliers off
jag_flat <- jag$clone(deep = TRUE) # necessary for r6, no overwriting!
jag_flat$multipliers[] <- 1L
r3 <- jag_flat$eval_pp(par = ssfit$par, n_jump = jag$max_dist - 1, env_type = "1o")
message(paste("PP at SS limit:", r3$ll_total,
              "| gap vs SS:", round(r3$ll_total - r1$ll_total, 4)))

