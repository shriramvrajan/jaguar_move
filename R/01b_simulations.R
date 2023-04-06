source("R/00_basics.R")

env1 <- gen_landscape(size = 100, s = 0.2, r = 30)
# env2 <- gen_landscape(size = 100, s = 0.033, r = 20)

x0 <- sample(100, 1); y0 <- sample(100, 1)
# x0 <- 50; y0 <- 50

path <- jag_path(x0, y0, 300, par = -1, neighb = 6)

plot_path(path, cut = 100)

### Ideas
# - Batch simulations to try and recover env AC from path AC
# - Home range weighting in path generation
# - Two turn angle distributions?

# Home range idea captures mechanism of other jaguars
# Coefficient could go to 0 in this case (driven by env AC)

# plan 1. compare AKDE ranges vs env suitability for real jaguars
# plan 2. simultaneously fit both env and home range parameters

# AKDE: pro, easy, con, maybe not logical

# How to give primacy to measured environmental variables?

# Preference not necessarily the same thing as suitability 

# Test cases.
# 1. XY GAM with threshold of occurrence
# 2. Two env layers and try to recover parameter values
# 3. What happens if you always start in a patch vs 'desert'

# Cardille 
# Jeff - modeling spatial challenges
# Melanie - behavioral ecology