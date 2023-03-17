## Basic functions and libraries

# Libraries
library(raster) # update to terra at some point
library(tidyverse)
library(data.table)
library(ctmm)
library(finch)
library(amt) 
library(mixtools)
library(gstat)

# Functions ====================================================================

# General ----------------------------------------------------------------------

# Save raster x as filename fn under ./data/
save_ras <- function(x, fn) {
    # Save raster often because it can get corrupted while working
    writeRaster(x, paste0("data/", fn), format = "raster",
                overwrite = TRUE)
}

# Output message m in logfile f
msg <- function(m, f = "data/output/run_log.txt") {
    m <- paste(format(Sys.time(), "%d.%m.%y  %R"), m)
    cat(m, file = f, append = TRUE, sep = "\n")
}

# Load files from a list if they exist in a directory dir
load_if_exists <- function(files, dir) {
    out <- lapply(files, function(f) {
        if (file.exists(paste0(dir, "/", f))) {
            readRDS(paste0(dir, "/", f))
        } else {
            return(NA)
        }
    })
}

# Simulation -------------------------------------------------------------------

# Outputs a matrix of cell numbers corresponding to raster (r, rdf)
# based on a central cell (i) and a buffer size around that cell (sz)
make_nbhd <- function(r = brazil_ras, rdf = brdf, i, sz) {
  # if x is the center square and o are neighbors, and e.g. sz = 2
  # (2*sz + 1)^2 represents the following:
  # o o o o o
  # o o o o o
  # o o x o o
  # o o o o o
  # o o o o o

  mat <- matrix(0, nrow = length(i), ncol = (2 * sz + 1)^2)

  # values to add to central cell's row/col to get neighborhood cells' row/col
  ind1 <- t(rep(-sz:sz, each = 2 * sz + 1))
  ind2 <- t(rep(-sz:sz, 2 * sz + 1))
  for (j in seq_len(length(ind1))) {
    mat[, j] <- cellFromRowCol(
      r, rdf$row[i] + ind1[j],
      rdf$col[i] + ind2[j]
    )
  }
  return(mat)
}

# Generates a random field with a given correlation structure
# b: beta parameter for gstat, s: sill, r: range, n: nugget
gen_landscape <- function(size = 100, b = 1, s = 0.03, r = 10, n = 0) {   
    xy <- expand.grid(1:size, 1:size); names(xy) <- c("x", "y")
    model <- gstat(formula = z ~ 1, locations = ~x + y, dummy = T, beta = b, 
                model = vgm(psill = s, range = r, model = "Exp",
                            nugget = n), nmax = 20)
    out <- predict(model, newdata = xy, nsim = 1)
    if (any(out < 0)) out[out < 0] <- 0
    gridded(out) <- ~x + y; out <- raster(out)
    raster::plot(out)

    outdf <- as.data.frame(out)
    outdf <- cbind(outdf, rowColFromCell(out, seq_len(nrow(outdf))))
    return(list(raster = out, df = outdf))
}

# Generate a jaguar path of n steps starting from (x1, y1) with environmental
# preference parameters par[] and search neighborhood size neighb
jag_path <- function(x1, y1, nstep, par = c(1, 1), neighb = 5) {
    if (!(x1 %in% 1:100) | !(y1 %in% 1:100)) {
        print("Jaguar out of bounds")
        return(NULL)
    }
    path <- matrix(NA, nrow = nstep, ncol = 2)
    x <- x1; y <- y1; path[1, ] <- c(x, y)
    for (i in 2:nstep) {
        pos <- path[i - 1, ]
        nbhd <- make_nbhd(r = env1[[1]], rdf = env1[[2]], sz = neighb,
                          i = cellFromRowCol(env1[[1]], pos[1], pos[2]))
        attract <- env1[[1]][nbhd] * par[1] + env2[[1]][nbhd] * par[2]
        if (any(is.na(attract))) attract[is.na(attract)] <- 0
        if (any(attract < 0)) attract[attract < 0] <- 0
        attract <- attract / sum(attract)
        step <- sample(seq_len(length(attract)), 1, prob = attract)
        path[i, ] <- rowColFromCell(env1[[1]], nbhd[step])
    }
    path <- as.data.frame(path)
    names(path) <- c("x", "y")
    return(path)
}

vgram <- function(path, cut = 100) {
    var <- sapply(1:cut, function(t) {
        p1 <- path[1:(100 - t),]
        p2 <- path[(t + 1):100,]

        out <- sqrt((p1$x - p2$x)^2 + (p1$y - p2$y)^2)
        return(mean(out))
    })
    return(var)
}

# Plot landscape r with jaguar path and vgram
plot_path <- function(path) {
    par(mfrow = c(1, 3))
    raster::plot(env1[[1]])
    points(path, col = "red", pch = 19, cex = 0.5)
    lines(path, col = "red")
    raster::plot(env2[[1]])
    points(path, col = "red", pch = 19, cex = 0.5)
    lines(path, col = "red")
    plot(vgram(path), type = "l", xlab = "Time lag", ylab = "Variance")
}



# Data =========================================================================

# WGS84 projection
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# Jaguar movement data, ID numbers, and metadata
jag_move <- readRDS("data/jag_data_BR.RDS")
jag_id <- readRDS("data/jag_list.RDS"); njag <- length(jag_id)
jag_meta <- data.table(read.csv("data/input/jaguars/jaguar_metadata.csv"))

# RasterStack of environmental variables 
# see 01_generate_data.R for details
brazil_ras <- stack("data/env_layers.grd")
# RasterStack of environmental variables, but as a data frame ('brdf')
load("data/env_layers.RData")

msg("Loaded data")

# ==============================================================================