library(sp)
library(sf)
library(terra)
source("./R/utils.R")

load("./clean_data/jaguar_trap_mats_js.rda")
load("./clean_data/grid_objs_data.rda")

## handle the environmental variables 
grid_pts_orig <- st_coordinates(grid_objs$grid_sf) * grid_objs$scaling
#steady_env <- grep(list.files("./clean_data/rasters", full.names = T), 
#                pattern='wlprai', invert=TRUE, value=TRUE) |>
#            lapply(terra::rast) |>
#            lapply(terra::extract, grid_pts_orig) |>
#            Reduce(f = cbind)
envX <- list.files("./clean_data/rasters", pattern = "wlprai", full.names = TRUE) |>
    lapply(terra::rast) |>
    lapply(terra::extract, grid_pts_orig) #|>
    #lapply(function(w,ww){
    #    cbind(w,ww)
    #}, steady_env)

good_grids <- envX |> Reduce(f = cbind) |> rowallgood()|> which()

envXlist <- envX |> lapply(`[` , good_grids, )

envX <- array(0, dim = c(length(envXlist), dim(as.matrix(envXlist[[1]]))))
means <- colMeans(as.matrix(envXlist[[1]]))
sds <- apply(as.matrix(envXlist[[1]]), 2, sd)
for(i in 1:length(envXlist)){
    envX[i, , ] <- normalizing( as.matrix( envXlist[[i]] ), means, sds)
}

## handle grid pts
grid_pts <- st_coordinates(grid_objs$grid_sf)[good_grids,]


## get dimension informations
capdim <- dim( jaguar_trap_mats$cap_mat )
M <- capdim[1]
n_trap <- capdim[2]
Kmax <- capdim[3]
T <- capdim[4]


JS_stan_data <- list(
    # dimensions
    M = M,
    n_trap = n_trap,
    Kmax = Kmax,
    T = T,

    # detection history
    y = jaguar_trap_mats$cap_mat,
    deploy = jaguar_trap_mats$trap_mat,

    # spatial locations
    X = st_coordinates(grid_objs$trap_sf), # trap locations
    n_grid = length(good_grids), # number of grid points
    grid_pts = grid_pts, # grid points

    # environmental variables
    n_env = dim(envX)[3], # number of environmental variables
    envX = envX # environmental variables
)

save(JS_stan_data, file = "./clean_data/js_stan_data.rda")
