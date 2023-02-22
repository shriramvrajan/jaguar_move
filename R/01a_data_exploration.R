# library(finch)
# library(OpenStreetMap)
# library(ggmap)

# ggmap::register_google(key = "AIzaSyC4Ewh1nlrRtvivd2qg2yXbr3lALhZv8a0")

### Data =======================================================================

## RasterStack of K's environmental variables
env <- rast("data/env_layers.grd")
## RasterStack of environmental variables, but as a data frame ('brdf')
# load("data/brazil.RData")

# Camtrap stuff
cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
gbif <- dwca_read("data/gbif")

print("Loaded data")

### Functions ==================================================================

plot_env <- function(layer, bounds, path, ras = env) {
    plot(crop(ras[[layer]], bounds), main = names(ras)[layer])
    lines(as.data.frame(path), col = rgb(1, 0, 0, 0.3), pch = 4)
}

# Map path of jaguar i
map_jag <- function(i, type = 2) {
    # type: 1 is satellite map using ggmap, 2 is env layers
    moves <- jag_move[ID == as.numeric(i)]

    # Get bounding box of all GPS measurements
    bbox <- raster::extent(data.frame(x = moves$longitude, y = moves$latitude))

    # Get total tracking period
    dates <- as.Date(sapply(moves$timestamp, function(dt) {
        strsplit(as.character(dt), " ")[[1]][1]
    }), format = "%m/%d/%y")
    period <- difftime(dates[length(dates)], dates[1])

    names(moves)[3:4] <- c("x", "y")
    path <- sp::SpatialPoints(coords = moves[, 3:4], sp::CRS("+init=epsg:4326"))
    # path2 <- as.data.frame(sp::spTransform(path, OpenStreetMap::osm()))

    if (type == 1) {
        bboxgg <- bbox[c(1, 3, 2, 4)]
        names(bboxgg) <- c("left", "bottom", "right", "top")
        map0 <- ggmap::get_map(location = bboxgg, maptype = "satellite")
        ggmap(map0) +
            geom_point(aes(x = x, y = y),
                data = as.data.frame(path),
                color = "red"
            )
    } else if (type == 2) {
        par(mfrow = c(2, 3))
        for (i in 1:6) {
            plot_env(i, bbox, path)
        }
    }
}

### Analyses ===================================================================

jag_move <- jag_move[, !c("study_name", "country")]
jag_moveid <- unique(jag_move$ID)




##### require('ctmm') ----------------------------------------------------------

jj <- jag_move[ID == as.numeric(jag_id[3])]
jj <- as.telemetry(jj, timeformat = "auto")
vjj <- variogram(jj); plot(vjj)
guess <- ctmm.guess(jj, interactive = F)
jjg <- ctmm.fit(jj, guess)
kde <- akde(jj, jjg)
plot(kde); plot(jj, add = T)

##### require('amt') -----------------------------------------------------------
# jag_tracks <- lapply(jag_moveid, function(id) {
#     id <- as.numeric(id)
#     path <- jag_move[ID == id]

#     path <- track(x = path$longitude, y = path$latitude,
#             timestamp = path$timestamp, id = path$ID)
# })
# jag <- jag_meta[, c("ID", "Sex", "Estimated.Age", "Weight")]
# jag <- jag[jag$ID %in% jag_moveid]
# names(jag) <- c("id", "sex", "age", "wt")
# jag$age <- as.numeric(jag$age)
# jag_steps <- lapply(jag_tracks, steps)



## Movebank stuff ==============================================================

# https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# https://cran.r-project.org/web/packages/move/move.pdf
# http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# searchMovebankStudies(x = "jaguar", login = login)
# bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
#                               login = login0)



## What are they actually fitting in  in state space models?
##     prob of being in a state or transition probabilities?
##     we need to know that we are doing it better
## Probability transition matrix model - consider all possible steps
##                                     - nonlinear routes
## Tr matrix for state space 
## True state 
## 
