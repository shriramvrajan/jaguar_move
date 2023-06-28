## Based on work by Kasia Johnson (2021) =======================================
## Generate data frames for jaguar movement analyses
## !IMPORTANT!: Run this through master.R

library(gfcanalysis)
library(sf)
library(WorldClimTiles)
library(raster) 
# All data generation was done using package 'raster', but it is recommended
# to use 'terra' for all subsequent analyses. 

### Following code chunks only need to be run once, and data saved.

### Producing environmental data for Brazil ====================================
# brazil_ras <- stack("data/brazil_ras_4.grd")
# merged <- readRDS("data/input/merged.RDS")
# m1 <- stack("data/input/bioclim_masked.grd")

## Finding and downloading the BIOCLIM tiles that make up Brazil ---------------
bra      <- getData("GADM", country = "BRA", level = 0)
tilename <- tile_name(bra)
tile     <- tile_get(tilename, "bio") 
# Makes list with all the data for each of the tiles for Brazil
merged   <- tile_merge(tile) # takes a while
saveRDS(merged, "data/merged.RDS")

## Crop to required range: need shapefile with boundaries of interest ----------
# Making sure projection of outline is same as projection of bioclim data
proj4string(spdf) == proj4string(brazil_ras)
Masking
msg("Masking bioclim tiles to Brazil's borders...")
m1 <- mask(merged, spdf) # takes a while 
m1 <- brick("data/input/bioclim_masked.grd")
# save_ras(m1, "input/bioclim_masked.grd")
# Cropping (using created mask to crop to the extent of Brazil)
msg("Cropping bioclim tiles to Brazil's borders")
brazil_ras <- crop(m1, extent(bra_shp)) # takes a while
msg("Saving initial brazil_ras.grd")
save_ras(brazil_ras, "input/brazil_ras.grd")

## Human footprint -------------------------------------------------------------
footprint             <- raster("data/input/wildareas-v3-2009-human-footprint.tif") 
footprint_rast        <- projectRaster(footprint, brazil_ras) 
# brazil_ras <- stack("data/input/brazil_ras.grd")
footprint_rast <- raster("data/input/footprint_raster.grd")
msg("Adding footprint layer")
brazil_ras <- addLayer(brazil_ras, footprint_rast) 
save_ras(brazil_ras, "input/brazil_ras.grd")

## Elevation and slope ---------------------------------------------------------
msg("Adding elevation and slope")
dem_ras <- raster("data/input/srtm30_dem.grd")
dem_ras <- crop(dem_ras, brazil_ras)
slope_ras <- terrain(dem_ras, opt = "slope")
brazil_ras <- addLayer(brazil_ras, dem_ras)
brazil_ras <- addLayer(brazil_ras, slope_ras)
save_ras(brazil_ras, "input/brazil_ras.grd")

# ## Forest cover --------------------------------------------------------------
brazil_ras <- brick("data/input/brazil_ras.grd")
# tiles <- calc_gfc_tiles(spdf) # tiles of forest cover that cover Brazil
# download_tiles(tiles,
#     output_folder = paste0(getwd(), "/data"),
#     images = "treecover2000",
#     dataset = "GFC-2019-v1.7"
# )
# Might have to rerun a couple times to get everything, connection times out
data_list <- list.files(path = "data/input", all.files = TRUE, full.names = TRUE)
hansen_list <- data_list[grep("Hansen", data_list)]
hansen <- lapply(hansen_list, raster) 
# transforming .tif into raster
hansen <- hansen[c(-10, -16, -17, -21)] 
# have to remove #s (faulty tif files): 10, 16, 17, 21
msg("Merging forest cover tiles")
forest_cover <- do.call(merge, hansen) # takes a long time
msg("Reprojecting forest cover")
forest_rast <- projectRaster(forest_cover, brazil_ras, method = "ngb")
save_ras(forest_rast, "input/forest_cover.grd")
forest_rast <- raster("data/input/forest_cover.grd")
msg("Resampling forest cover layer")
forest_rast <- resample(forest_rast, brazil_ras)
msg("Adding forest cover layer")
brazil_ras <- addLayer(brazil_ras, forest_rast)
msg("Saving brazil_ras with forest cover")
save_ras(brazil_ras, "input/brazil_ras.grd")

# ## Modifications by Shriram Varadarajan (2022)================================

# ## Adding distance to water & road
# brazil_ras <- stack("data/input/brazil_ras.grd")
# brazil_ras <- subset(brazil_ras, c(20:23))
# names(brazil_ras) <- c("footprint", "elevation", "slope", "forestcover")
# save_ras(brazil_ras, "input/brazil_ras.grd")
brazil_ras <- brick("data/input/brazil_ras.grd")

## Distance from water raster layer
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.71c6r
msg("Adding distance to water layer")
water <- raster("data/input/streams/distance2water_30arcsec.tif")
water <- resample(water, brazil_ras, method = "bilinear")
brazil_ras <- addLayer(brazil_ras, water)

## Get road layer, produce distanceFrom raster layer
## Distance calculated in QGIS using 'proximity' function
dist <- raster("data/input/roads/dist_to_road.tif")
dist <- projectRaster(dist, brazil_ras)
msg("Adding distance from road layer")
brazil_ras <- addLayer(brazil_ras, dist)
msg("Saving input raster")
names(brazil_ras)[5:6] <- c("distwater", "distroad")
save_ras(brazil_ras, "env_layers.grd")

msg("Transforming raster to data frame (brdf)")
brdf <- as.data.frame(brazil_ras)
msg("Adding XY data to brdf")
brdf <- cbind(brdf, rowColFromCell(brazil_ras, 1:ncell(brazil_ras))) 
msg("Adding indices to brdf")
brdf$index <- seq_len(nrow(brdf))
msg("Saving brdf")
save(brdf, file = "data/env_layers.RData")

## Adding biome ================================================================

library(terra)
# Eventually will shift to terra for previous steps as well

biome <- vect("data/input/Brazil_biomes/Brazil_biomes.shp")
jags <- readRDS("data/jag_data_BR.RDS")
jags <- do.call(rbind, by(jags, jags$ID, function(x) x[1, ]))
jags$biome <- terra::extract(biome, jags[, 3:4])$name

jag_meta <- data.table(read.csv("data/input/jaguars/jaguar_metadata.csv"))
jag_meta <- merge(jag_meta, jags[, c("ID", "biome")], by = "ID")
jag_meta$biome[42] <- "Mata Atlântica" # fix for one jaguar, Argentina border
write.csv(jag_meta, "data/input/jaguars/jaguar_metadata.csv", row.names = FALSE)

msg("Environmental layers saved, end of generate_data.R.")