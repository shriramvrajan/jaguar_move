source("R/00_functions.R")

### Data =======================================================================

jag_move <- jag_move[, !c("study_name", "country")]
jag_moveid <- unique(jag_move$ID)
jag <- jag_meta[, c("ID", "Sex", "Estimated.Age", "Weight")]
jag <- jag[jag$ID %in% jag_moveid]
names(jag) <- c("id", "sex", "age", "wt")
jag$age <- as.numeric(jag$age)

env <- rast("data/env_layers.grd")

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

### Analyses ==================================================================
n <- 80

tr <- jag_track(jag_id[n])
st <- steps(tr)

summarize_speed(tr)
summarize_sl(tr)

diff(sampling_period(tr))

map_jag(jag_id[n])
map_homerange(jag_id[n])

## Movebank stuff ==============================================================

# https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# https://cran.r-project.org/web/packages/move/move.pdf
# http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# searchMovebankStudies(x = "jaguar", login = login)
# bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
#                               login = login0)


