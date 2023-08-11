library(sp)
library(terra)
library(sf)
# We are at UTM-17N, 
## below for projection from WGS84, this is essential as distance is a thing here
crs_lat_long <- CRS("+proj=longlat +datum=WGS84")

if(!file.exists("./clean_data/shoreline.tif")){
  lclu <- rast("./data/Guanecaste/shore_line/shoreline.tif") |> 
    terra::project(y = "+proj=utm +zone=17")
  terra::writeRaster(lclu,"./clean_data/shoreline.tif")
}

camera_sites <- list.files("./data/Guanecaste", 
pattern = "trap_loc", full.names = T) |>
  lapply(read.csv) |>
  lapply(function(w, crs_lat_long){
    SpatialPointsDataFrame(data.frame(x = w$Trap.Station.Longitude, 
                                y = w$Trap.Station.Latitude),
                                data = w, 
                                proj4string = crs_lat_long)
  }, crs_lat_long) |>
  lapply(spTransform, CRS("+proj=utm +zone=17"))
# only one year
camera_sites <- camera_sites
save(camera_sites, file = "./clean_data/Guanecaste/CT_loc.rda")

jaguar <- read.csv(".//data/Guanecaste/jaguar_capture.csv") 
jaguar <- SpatialPointsDataFrame(data.frame(x = jaguar$Trap.Station.Longitude, 
                                            y = jaguar$Trap.Station.Latitude),
                                 data = jaguar,
                                 proj4string = crs_lat_long) |>
  spTransform(CRS("+proj=utm +zone=17")) 
  
save(jaguar,file = "./clean_data/Guanecaste/jaguar.rda")
