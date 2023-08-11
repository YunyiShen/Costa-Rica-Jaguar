library(sp)
library(terra)

# hard code all the rasters
raw_raster_loc <- "./data/Guanecaste/SpatialCovariates"
raster_target <- "./clean_data/Guanecaste/rasters"
raw_raster_files <- c(
    #"dist_prkbndry/hdr.adf",
    #"dist_station/hdr.adf",
    #"dist_trail/hdr.adf",
    #"riverdensity2/hdr.adf",
    "distance_to_shore.tif"
)

raster_files_target <- c(
    #"dist_prkbndry.tif",
    #"dist_station.tif",
    #"dist_trail.tif",
    #"river_density.tif",
    "distance_to_shore.tif"
)

raw_raster_files <- file.path(raw_raster_loc, raw_raster_files)
raster_files_target <- file.path(raster_target, raster_files_target)

for(i in 1:length(raw_raster_files)){
    if(!file.exists(raster_files_target[i])){
        r <- terra::rast(raw_raster_files[i]) |> 
                terra::project(y = "+proj=utm +zone=17")
        terra::writeRaster(r,raster_files_target[i])
    }
}





