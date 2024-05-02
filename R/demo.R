source("R/00_functions.R")

envsize <- 40
s1      <- 5
r1      <- 5
sim_interval <- 0

# Landscape
set.seed(7)
env01 <- gen_landscape(envsize, s1, r1)
terra::plot(env01[[1]])

# Attraction function
par0         <- c(2, 0, -0.1) 
# plot_the_curve(par0)

# trad model: sample all cells
ext <- c(14, 25, 14, 25)
x0 <- ceiling(envsize / 2) ## Starting xy
y0 <- ceiling(envsize / 2)
terra::plot(env01[[1]], ext = ext)
points(x0, y0, col = "red", pch = 19, cex = 2)

# new model: one step at a time
par(mfrow = c(1, 1))
step1 <- expand.grid(x = (x0 - 1):(x0 + 1), y = (y0 - 1):(y0 + 1))
points(step1, col = "blue", pch = 19, cex = 2)

step2 <- expand.grid(x = (x0 - 2):(x0 + 2), y = (y0 - 2):(y0 + 2))
points(step2, col = "blue", pch = 19, cex = 2)

### shortsightedness -----------------------------------------------------------
n_step <- 200
step_size    <- 1
set.seed(16)
p1 <- jag_path(x0, y0, n_step, par = par0, neighb = step_size)
plot_path(p1)
points(x0, y0, col = "red", pch = 19, cex = 2)
points(step1, col = "blue", pch = 19, cex = 2)

n_step <- 100
step_size    <- 3
set.seed(10)
p2 <- jag_path(x0, y0, n_step, par = par0, neighb = step_size)
plot_path(p2)
points(x0, y0, col = "red", pch = 19, cex = 2)


## home range ------------------------------------------------------------------

# homerange_41.png
env <- brazil_ras

plot_xyt(make_track0(41))

## 3d plot
tr <- make_track0(20)
plot_xyt(tr)

## multiple dispersal kernels

# homerange_41.png again - noise and signal
