load_source <- TRUE

if (load_source) {
    source("R/03_results.R")
}

plot_indiv <- TRUE
 
### Data =======================================================================

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

### Analyses ===================================================================

# Testing holdout sets
i = 2
message(paste0("Jaguar #: ", i))
id <- as.numeric(jag_id[i])

jag_traject <- jag_move[ID == id, 3:4]
holdout_frac <- 0.3 
hold <- seq_len(ceiling(nrow(jag_traject) * holdout_frac))
jag_traject <- jag_traject[-hold, ]


# Plotting for individual jaguar -----------------------------------------------
if (plot_indiv) {
    n <- 20
    map_track(n)

    j <- add_track_metadata(n)
    hist(j$ta)

    # ctmm + plotly stuff for 3D plot ------------------------------------------
    tel <- as.telemetry(jag_move[ID == n], timeformat = "auto")
    # lunarize(tel, 1)
    fig <- plot_ly(x = tel$x, y = tel$y, z = tel$t, type = "scatter3d", 
                   mode = "lines", line = list(width = 1)) 
    # km <- kmeans(tel[, 4:6], centers = 1)
    # tel$clus <- km$cluster
    # fig <- plot_ly(x = tel$x, y = tel$y, z = tel$t, type = "scatter3d", 
    #                mode = "lines", line = list(width = 1, color = tel$clus)) %>% 
    #                add_trace(x = km$centers[, 2], y = km$centers[, 3], 
    #                          z = km$centers[, 1], type = "scatter3d", 
    #                          mode = "markers", marker = list(size = 5),
    #                          line = list(color = "white"))
    
    # kmeans stuff -------------------------------------------------------------
    wss <- vector()
    for (k in seq_len(10)) {
        km <- kmeans(tel[, 4:6], centers = k)
        wss[k] <- km$tot.withinss
    }
    plot(wss)

}

## Movebank stuff ==============================================================

# https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# https://cran.r-project.org/web/packages/move/move.pdf
# http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# searchMovebankStudies(x = "jaguar", login = login)
# bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
#                               login = login0)