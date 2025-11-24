library(devtools)
library(WorldClimTiles)
library(rgdal)
library(maptools)
library(geosphere)
library(sp)
library(adehabitatHR)
library(gfcanalysis)
library(utils)
library(ncdf4)

### Producing environmental data for Brazil ----

# Finding and downloading the BIOCLIM tiles that make up Brazil:
bra <- getData("GADM", country = "BRA", level = 0)
tilename <- tile_name(bra)
tile <- tile_get(tilename, "bio") # Makes list with all the data for each of the tiles for Brazil
merged <- tile_merge(tile) # takes a while

# Crop to required range: need a shapefile that denotes the boundaries of interest
data("wrld_simpl")
spdf <- subset(wrld_simpl, NAME == "Brazil") # spatial polygon of Brazil
spdf <- spTransform(spdf, CRSobj = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
# # making sure projection of Brazil outline is the same as the projection of bioclim data
proj4string(spdf) == proj4string(merged)
plot(spdf)

# # Masking:
m1 <- mask(merged, spdf) # takes a while

# Cropping (using created mask to crop to the extent of Brazil):
brazil_ras <- crop(m1, extent(spdf)) # takes a while
footprint <- raster("wildareas-v3-2009-human-footprint.tif") # Human footprint
footprint_rast <- projectRaster(footprint, brazil_ras) # projecting to Brazil area
ncell(footprint_rast)==ncell(brazil_ras)

saveRDS(m1, "data/m1.RDS")
saveRDS(brazil_ras, "data/brazil_ras.RDS")

merged <- readRDS("merged.RDS")
m1 <- readRDS("m1.RDS")
brazil_ras <- readRDS("brazil_ras.RDS")
footprint_rast <- readRDS("footprint_rast.RDS")




## Adding additional variables not included in bioclim ---

# Human footprint index:

brazil_ras <- addLayer(brazil_ras, footprint_rast) # adding to brazil_ras RasterBrick

# Forest cover:

tiles <- calc_gfc_tiles(spdf) # tiles of forest cover that cover Brazil
# plot(tiles)
# plot(spdf, add=T)
download_tiles(tiles,
    output_folder = getwd(),
    images = "treecover2000",
    dataset = "GFC-2019-v1.7"
)


setwd("~/Desktop/Scripts/Hansen")
hansen_list <- list.files(all.files = T, full.names = T) 
# also returned two empty values ("./.)
hansen_list <- hansen_list[c(-1, -2)]
hansen <- lapply(hansen_list, raster) 
# transforming .tif into raster
hansen <- hansen[c(-10, -16, -17, -21)] 
# have to remove #s (faulty tif files): 10, 16, 17, 21



Sys.time()
forest_cover <- do.call(merge, hansen) # takes a long time
Sys.time()
beep(5)
saveRDS(forest_cover, file = "forest_cover.rds")

forest_rast <- projectRaster(forest_cover, brazil_ras, method = "ngb")
brazil_ras <- addLayer(brazil_ras, forest_rast)





# Elevation and Slope:

setwd("~/Desktop/Scripts")
DEM <- raster("srtm30plus_v11_land.nc")
DEM_ras <- projectRaster(DEM, brazil_ras, method = "ngb")
slope_ras <- terrain(DEM_ras, opt = "slope")
brazil_ras <- addLayer(brazil_ras, DEM_ras) # adding to brazil_ras RasterBrick
brazil_ras <- addLayer(brazil_ras, slope_ras)

# Dropping layers not needed to save space

brazil_ras <- dropLayer(brazil_ras, c(2:11, 13:19))





# Convert to data frame & save
brdf <- as.data.frame(brazil_ras)
brdf <- cbind(brdf, rowColFromCell(brazil_ras, 1:ncell(brazil_ras)))
brdf$index <- 1:nrow(brdf)





# 	save(brdf, file = "brazil.rds")
saveRDS(brdf, "brazil.rds")
saveRDS(brazil_ras, "brazil_raster.rds")
