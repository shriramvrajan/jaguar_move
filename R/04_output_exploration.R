plot_res <- FALSE
plot_indiv <- TRUE
load_source <- TRUE

if (load_source) {
    source("R/03_results.R")
}

### Data =======================================================================

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

### Analyses ===================================================================


# Plotting from results
if (plot_res) {
    # vars: 1 footprint, 2 elev, 3 slope, 4 forestcover, 5 distwater, 
    #       6 distroad, 7 homerange
    barplot(tapply(res$p1, res$bio, function(x) mean(x, na.rm = TRUE)),
            col = "#a37373", border = NA)

    evars <- c("footprint", "elev", "slope", "forest", "distwater", "distroad",
            "homerange")

    res <- as.data.frame(res) # still not used to plotting with tibbles

    

    par(mfrow = c(2, 4))

    for (i in seq_len(length(evars))) {
        hist(res[, paste0("p", i)], 20, col = "#a37373", border = NA,
             xlab = "Parameter value", main = evars[i])
        abline(v = 0, lwd = 2, col = "blue")
        print(evars[i])
        print(summary(res[, paste0("p", i)]))
    }
}


# Plotting for individual jaguar
if (plot_indiv) {
    n <- which(jag_id == 20)

    tr <- jag_track(jag_id[n])
    st <- steps(tr)

    print(summarize_speed(tr))
    print(summarize_sl(tr))

    print(diff(sampling_period(tr)))

    # ggplot(aes(x = t_, y = y_), data = tr)
    plot_ly(x = x_, y = y_, z = t_, type = "scatter3d", mode = "markers", color = temp)
    # print(map_homerange(jag_id[n]))
}

## Movebank stuff ==============================================================

# https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# https://cran.r-project.org/web/packages/move/move.pdf
# http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# searchMovebankStudies(x = "jaguar", login = login)
# bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
#                               login = login0)