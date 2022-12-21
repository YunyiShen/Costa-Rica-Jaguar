library(sp)
library(terra)
library(sf)
# We are at UTM-17N, 
## below for projection from WGS84, this is essential as distance is a thing here
crs_lat_long <- CRS("+proj=longlat +datum=WGS84")

if(!file.exists("./clean_data/lulc_2017.tif")){
  lclu <- rast("./data/mask/NASA_OSA/lulc_2017.tif") |> 
    terra::project(y = "+proj=utm +zone=17")
  terra::writeRaster(lclu,"./clean_data/lulc_2017.tif")
}

camera_sites <- list.files("./data/Camera_loc", pattern = "CostaCT", full.names = T) |>
  lapply(read.csv) |>
  lapply(function(w, crs_lat_long){
    SpatialPointsDataFrame(data.frame(x = w$gps_x, y = w$gps_y),
                                data = w, 
                                proj4string = crs_lat_long)
  }, crs_lat_long) |>
  lapply(spTransform, CRS("+proj=utm +zone=17")) 
save(camera_sites, file = "./clean_data/CT_loc.rda")

jaguar <- read.csv("./data/detections/JaguarEvents2015-2021_EDIT.csv") 
jaguar <- SpatialPointsDataFrame(data.frame(x = jaguar$Trap.Station.Longitude, 
                                            y = jaguar$Trap.Station.Latitude),
                                 data = jaguar,
                                 proj4string = crs_lat_long) |>
  spTransform(CRS("+proj=utm +zone=17")) 
  
save(jaguar,file = "./clean_data/jaguar.rda")
