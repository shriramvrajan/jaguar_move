source("R/00_basics.R")

### Functions ==================================================================

### Generating two random fields with different correlation structures =========

env1 <- gen_landscape(size = 100, s = 0.02, r = 15)
env2 <- gen_landscape(size = 100, s = 0.003, r = 20)

# x0 <- sample(100, 1); y0 <- sample(100, 1)
x0 <- 50; y0 <- 50

path <- jag_path(x0, y0, 100, par = c(1, 0), neighb = 6)

plot_path(path)

### Next steps: ================================================================