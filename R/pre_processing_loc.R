library(sp)
library(terra)
library(sf)
config <- jsonlite::fromJSON("config.json")
# We are at UTM-17N, 
## below for projection from WGS84, this is essential as distance is a thing here
crs_lat_long <- CRS("+proj=longlat +datum=WGS84")

if(!file.exists("./clean_data/lulc_2017.tif")){
  lclu <- rast("./data/mask/NASA_OSA/lulc_2017.tif") |> 
    terra::project(y = "+proj=utm +zone=17")
  terra::writeRaster(lclu,"./clean_data/lulc_2017.tif")
}

camera_sites <- list.files("./data/Camera_loc", pattern = "CostaCT", full.names = T) |>
  lapply(read.csv) 

for(i in 1:length(camera_sites)){
  names(camera_sites[[i]]) <- c("Station", "Camera","gps_y", "gps_x", "Setup_date", "Retrieval_date") 
}
camera_sites <- camera_sites  |>
  lapply(function(w, crs_lat_long){
    st_as_sf(
      w, coords = c("gps_x", "gps_y"), crs = crs_lat_long
    )
  }, crs_lat_long) |>
  lapply(st_transform, CRS("+proj=utm +zone=17")) 

for(i in 1:length(camera_sites)){
  camera_sites[[i]]$x <- st_coordinates(camera_sites[[i]])[,1]
  camera_sites[[i]]$y <- st_coordinates(camera_sites[[i]])[,2]
}

save(camera_sites, file = "./clean_data/CT_loc.rda")

jaguar <- read.csv(config$detection) 
jaguar <- st_as_sf(
  jaguar, coords = c("Trap.Station.Longitude", "Trap.Station.Latitude"), crs = crs_lat_long
  ) |>
  st_transform(CRS("+proj=utm +zone=17")) 
xy_jaguars <- st_coordinates(jaguar)
jaguar$x <- xy_jaguars[,1]
jaguar$y <- xy_jaguars[,2]


save(jaguar,file = "./clean_data/jaguar.rda")
