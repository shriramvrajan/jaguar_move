## Basic functions and libraries

# Libraries
library(raster) # update to terra at some point
library(tidyverse)
library(data.table)
library(ctmm)
library(finch)
library(amt) 

# Functions ====================================================================

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