load_source <- TRUE

if (load_source) {
    source("R/03_results.R")
}

## Set up parallel processing ==================================================

ncore <- 2

library(doParallel)
library(foreach)
registerDoParallel(ncore)
message(paste0("number of workers: ", getDoParWorkers()))


### Analyses ===================================================================

# Testing holdout sets
# A lot of this code is also in 02_movement_model ; refactor later


jag <- seq_len(nrow(jag_id))

bad <- c(1, 2, 33, 42, 49)

# jag <- jag[1:8]

# fitholdout <- foreach(i = jag, .combine = "rbind") %dopar% {
fitholdout <- lapply(jag, function(i) {
    # browser()
    message(paste0("Jaguar #: ", i))
    
    if (i %in% bad) {
        return(cbind(NA, NA))
    }

    id <- as.numeric(jag_id[i])

    nbhd0 <- make_nbhd(i = seq_len(nrow(brdf)), sz = buffersize)
    jag_traject <- jag_move[ID == id, 3:4]
    holdout_frac <- 0.3
    hold <- seq_len(ceiling(nrow(jag_traject) * holdout_frac))
    jag_traject <- jag_traject[-hold, ]

    jag_traject_cells <- cellFromXY(brazil_ras, jag_traject)
    n_obs <- length(jag_traject_cells)
    # Calculating step distances; divide by cell size then take hypotenuse
    dist <- (jag_traject[-nrow(jag_traject), ] - jag_traject[-1, ]) /
        xres(brazil_ras)
    dist <- (rowSums(dist^2))^.5
    max_dist <- ceiling(max(dist) * 2)
    prep_model_objects(jag_traject_cells, max_dist,
        r = brazil_ras, rdf = brdf,
        nbhd0 = nbhd0
    )
    sim_steps <- 25

    # Normalizing desired environmental variables for extended neighborhood
    home <- rast(paste0("data/homeranges/homerange_", id, ".grd"))
    brdf$home <- as.vector(home)
    envdf <- brdf[, c(1:6, 10)]
    env <- envdf[nbhd_index, ]
    env <- sweep(env, 2, colMeans(env), "-") 
    env <- sweep(env, 2, apply(env, 2, sd), "/") 
    # Make indexing consistent with env
    row.names(env) <- seq_len(length(nbhd_index))
    
    print("done prep")

    par1 <- param$holdRWM[i, ]
    par2 <- param$holdtrad1[i, ]
    # sims <- paste0("data/output/", c("LL_holdRWH", "LL_holdtrad1"))
    # par1 <- readRDS(paste0(sims[1], "/par_out_", i, ".RDS"))
    # par2 <- readRDS(paste0(sims[2], "/par_out_", i, ".RDS"))
    objects <- list(env, nbhd, max_dist, n_obs, sim_steps, to_dest, obs)
    ll_1 <- log_likelihood(par1, objects)
    print("done ll1")

    track <- make_full_track(id)
    sl_emp <- as.vector(na.exclude(track$sl))
    ta_emp <- as.vector(na.exclude(track$ta))
    mk <- make_movement_kernel(sl_emp, ta_emp, n = 10000, max_dist = max_dist)
    objects <- list(env, max_dist, mk, obs)
    ll_2 <- log_likelihood0(par2, objects)
    print("done ll2")
    return(cbind(ll_1, ll_2))
})
fitholdout <- do.call(rbind, fitholdout)

saveRDS(fitholdout, "data/output/fitholdout.RDS")

hold <- as.data.frame(readRDS("data/output/fitholdout.RDS"))
plot(hold$ll_1, hold$ll_2)
abline(0, 1)

aic1 <- 2 * 7 - 2 * (-hold$ll_1)
aic2 <- 2 * 7 - 2 * (-hold$ll_2)

### Data =======================================================================

# Camtrap stuff
# cam_unit <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_UNIT.csv")
# cam_area <- read.csv("data/camtrap/dataset/AMZ_CAMTRAP_AREA.csv")

# GBIF records
# gbif <- read.csv("data/gbif_jaguars_brazil.csv")
# gbif <- dwca_read("data/gbif")

# Plotting for individual jaguar -----------------------------------------------
# if (plot_indiv) {
#     n <- 20
#     map_track(n)

#     j <- add_track_metadata(n)
#     hist(j$ta)

#     # ctmm + plotly stuff for 3D plot ------------------------------------------
#     tel <- as.telemetry(jag_move[ID == n], timeformat = "auto")
#     # lunarize(tel, 1)
#     fig <- plot_ly(
#         x = tel$x, y = tel$y, z = tel$t, type = "scatter3d",
#         mode = "lines", line = list(width = 1)
#     )
#     # km <- kmeans(tel[, 4:6], centers = 1)
#     # tel$clus <- km$cluster
#     # fig <- plot_ly(x = tel$x, y = tel$y, z = tel$t, type = "scatter3d",
#     #                mode = "lines", line = list(width = 1, color = tel$clus)) %>%
#     #                add_trace(x = km$centers[, 2], y = km$centers[, 3],
#     #                          z = km$centers[, 1], type = "scatter3d",
#     #                          mode = "markers", marker = list(size = 5),
#     #                          line = list(color = "white"))

#     # kmeans stuff -------------------------------------------------------------
#     wss <- vector()
#     for (k in seq_len(10)) {
#         km <- kmeans(tel[, 4:6], centers = k)
#         wss[k] <- km$tot.withinss
#     }
#     plot(wss)
# }

# ## Movebank stuff ==============================================================

# # https://cran.r-project.org/web/packages/move/vignettes/browseMovebank.html
# # https://cran.r-project.org/web/packages/move/move.pdf
# # http://www2.uaem.mx/r-mirror/web/packages/move/vignettes/move.pdf

# # login0 <- movebankLogin(username = "shriramv", password = "nY@tR6YT")

# # searchMovebankStudies(x = "jaguar", login = login)
# # bci_ocelot <- getMovebankData(study = "Ocelots on Barro Colorado Island, Panama",
# #                               login = login0)
