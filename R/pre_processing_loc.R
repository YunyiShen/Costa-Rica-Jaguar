library(sp)
library(terra)
# We are at UTM-17N, 
## below for projection from WGS84, this is essential as distance is a thing here
crs_lat_long <- CRS("+proj=longlat +datum=WGS84")
camera_sites <- list.files("./data", pattern = "CostaCT", full.names = T) |>
  lapply(read.csv) |>
  lapply(function(w, crs_lat_long){
    SpatialPointsDataFrame(data.frame(x = w$gps_x, y = w$gps_y),
                                data = w, 
                                proj4string = crs_lat_long)
  }, crs_lat_long)

jaguar <- read.csv("./data/JaguarEvents2015-2021_EDIT.csv") 
jaguar <- SpatialPointsDataFrame(data.frame(x = jaguar$Trap.Station.Longitude, 
                                            y = jaguar$Trap.Station.Latitude),
                                 data = jaguar,
                                 proj4string = crs_lat_long)
