library(sp)
library(sf)
library(terra)
library(jsonlite)
source("./R/utils.R")

load("./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")
lclu <- terra::rast("./clean_data/shoreline.tif")
config <- jsonlite::fromJSON("config.json")

scaling <- config$scaling
# setup a grid to approximate the marginalization over space
# smaller delta values --> better approximation
delta <- config$resolution
buffer <- config$buffer


trap_mat <- jaguar_trap_mats$trap_mat
traplocs <- jaguar_trap_mats$ids$trap_ids[,c("x","y")] / scaling # km's
K.jaguar <- rowSums(trap_mat) # row sums
cap_mat <- jaguar_trap_mats$cap_mat


#y3d <- SCR23darray(cap_mat, trap_mat)
#y2d <- apply(y3d, c(1, 2), sum)



x1_grid <- seq(min(traplocs$x) - buffer, 
               max(traplocs$x) + buffer, 
               by = delta)
x2_grid <- seq(min(traplocs$y) - buffer, 
               max(traplocs$y) + buffer, 
               by = delta)
grid_pts <- expand.grid(x = x1_grid, y = x2_grid) 
lu_grid <- terra::extract(y = grid_pts*scaling, x = lclu)
grid_pts <- grid_pts[!is.na(lu_grid$shoreline),]


# filter points based on distance to traps ---------------
grid_sf <- st_as_sf(grid_pts, coords = c("x", "y"))
trap_sf <- st_as_sf(traplocs, coords = c("x", "y"))
trap_buffer <- st_buffer(trap_sf, buffer)

# this horrorshow removes grid points outside trap buffers 
grid_sf <- grid_sf[unlist(lapply(st_intersects(grid_sf, trap_buffer), length)) > 0, ]

#plot(grid_sf, pch = 4, col = "red", cex = .3)
#plot(trap_sf, add = TRUE)

grid_objs <- list(grid_sf = grid_sf, trap_sf = trap_sf, trap_buffer = trap_buffer, traplocs = traplocs, scaling = scaling)
save(grid_objs, file = "./clean_data/Guanecaste/grid_objs_data.rda")
