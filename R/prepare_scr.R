# get src ready 
library(sp)
library(sf)
library(terra)
library(jsonlite)
source("./R/utils.R")
config <- fromJSON("./config.json")

load("./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")
load("./clean_data/Guanecaste/grid_objs_data.rda")

trap_mat <- jaguar_trap_mats$trap_mat
cap_mat <- jaguar_trap_mats$cap_mat

grid_pts_orig <- st_coordinates(grid_objs$grid_sf) * grid_objs$scaling

envX <- list.files("./clean_data/Guanecaste/rasters", pattern = config$env_key_words, full.names = TRUE) |>
    lapply(terra::rast) |>
    lapply(terra::extract, grid_pts_orig) |>
    Reduce(f = cbind) |>
    as.matrix()

good_grids <- envX  |> rowallgood()|> which()

envX <- envX[good_grids,] |>
    as.matrix()
means <- colMeans(as.matrix(envX))
sds <- apply(as.matrix(envX), 2, sd)


envX <- apply(envX,1,function(x) (x-means)/sds) |>
  as.matrix()


## handle grid pts
grid_pts <- st_coordinates(grid_objs$grid_sf)[good_grids,]


## get dimension informations
capdim <- dim( jaguar_trap_mats$cap_mat )
M <- capdim[1]
n_trap <- capdim[2]
Kmax <- capdim[3]

scr_stan_data <- list(
    # dimensions
    M = M,
    n_trap = n_trap,
    Kmax = Kmax,

    # detection history
    y = jaguar_trap_mats$cap_mat,
    deploy = jaguar_trap_mats$trap_mat,

    # spatial locations
    X = st_coordinates(grid_objs$trap_sf), # trap locations
    n_grid = length(good_grids), # number of grid points
    grid_pts = grid_pts, # grid points

    # environmental variables
    n_env = dim(envX)[2], # number of environmental variables
    envX = envX # environmental variables
)

save(scr_stan_data, file = "./clean_data/Guanecaste/scr_stan_data.rda")



