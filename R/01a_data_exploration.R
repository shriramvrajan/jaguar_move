source("R/00_functions.R")

### Data =======================================================================

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

env <- rast("data/env_layers.grd")

jag <- jag_meta[ID %in% jag_move$ID]

### Analyses ===================================================================
n <- 56

tr <- jag_track(jag_id[n])
st <- steps(tr)

print(summarize_speed(tr))
print(summarize_sl(tr))

print(diff(sampling_period(tr)))

print(map_jag(jag_id[n]))
print(map_homerange(jag_id[n]))

## Movebank stuff ==============================================================

# https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# https://cran.r-project.org/web/packages/move/move.pdf
# http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# searchMovebankStudies(x = "jaguar", login = login)
# bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
#                               login = login0)