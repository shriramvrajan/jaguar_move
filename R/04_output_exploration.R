load_source <- TRUE

if (load_source) {
    source("R/03_results.R")
}

plot_res <- FALSE
plot_indiv <- TRUE
 
### Data =======================================================================

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

### Analyses ===================================================================

# Plotting from results --------------------------------------------------------
if (plot_res) {
    # vars: 1 footprint, 2 elev, 3 slope, 4 forestcover, 5 distwater, 
    #       6 distroad, 7 homerange
    evars <- c("footprint", "elev", "slope", "forest", "distance to water", "distance to road",
            "unfamiliarity")

    i <- 3
    
    barplot(tapply(res[, paste0("p", i)], res$bio, function(x) mean(x, na.rm = TRUE)),
            col = "#a37373", border = NA, main = evars[i])



    res <- as.data.frame(res) # still not used to plotting with tibbles

    

    par(mfrow = c(2, 4))

    for (i in seq_len(length(evars))) {
        hist(res[, paste0("p", i)], 20, col = "#a37373", border = NA,
             xlab = "Parameter value", main = evars[i])
        abline(v = 0, lwd = 2, col = "black")
        print(evars[i])
        print(summary(res[, paste0("p", i)]))
    }
}


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
    fig
    
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