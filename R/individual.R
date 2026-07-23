source("R/functions.R")
source("R/classes.R")

test_warmstart <- FALSE

ss_res <- readRDS("data/output/emp_ss_1o.rds")
pp_res <- readRDS("data/output/emp_pp_2026-07-17.rds")
r1 <- results_set$new(r_ss = "data/output/emp_ss_1o.rds", 
        r_pp = "data/output/emp_pp_2026-07-17.rds", env_type = "1o")$res_table

id  <- 113
jag <- jaguar$new(id)
par0 <- ss_res[[which(jag_id$jag_id == id)]]$par
par1 <- pp_res[[which(jag_id$jag_id == id)]]$par

if (test_warmstart) {
    result    <- jag$fit_pp(par0, env_type = "1o")
    env_index <- head(seq_along(par0), -2) # -last 2
    par1 <- c(rnorm(length(env_index), 0, 0.1), tail(par0, 2))
    result2   <- jag$fit_pp(par1, env_type = "1o")
    saveRDS(list("warm" = result, "cold" = result2), "data/output/test.rds")
}

# result2   <- jag$fit_pp(par1, env_type = "1o", n_jump = 5, clamp = Inf)
# par2 <- c(2.56, -1.15, -.38, -14.2, 5.72, -7.5, -9.7, 0.77, -9.88)

ssfit <- jag$fit_ss(par0)
ppfit <- jag$fit_pp(par0)

r1 <- jag$eval_ss(par = ssfit$par, env_type = "1o")
message(paste("SS log(L):", r1$ll_total))

## PP at the grain the sweep chose
r2 <- jag$eval_pp(par = ppfit$par, n_jump = ppfit$n_jump)
message(paste("PP log(L):", r2$ll_total, "| n_jump =", ppfit$n_jump,
              "| best_r =", r2$best_r))

## nesting check: multipliers off, one hop spans the window
jag_flat <- jag$clone(deep = TRUE) # necessary for r6, no overwriting!
jag_flat$multipliers[] <- 1L
r3 <- jag_flat$eval_pp(par = ssfit$par, n_jump = jag$max_dist - 1, env_type = "1o")
message(paste("PP at SS limit:", r3$ll_total,
              "| gap vs SS:", round(r3$ll_total - r1$ll_total, 4)))