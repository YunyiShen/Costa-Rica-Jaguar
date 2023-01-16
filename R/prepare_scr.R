# get src ready 
library(sp)
library(sf)
library(terra)
source("./R/utils.R")

load("./clean_data/jaguar_trap_mats_scr.rda")
load("./clean_data/grid_objs_data.rda")

trap_mat <- jaguar_trap_mats$trap_mat
K.jaguar <- rowSums(trap_mat) # row sums
cap_mat <- jaguar_trap_mats$cap_mat


y3d <- SCR23darray(cap_mat, trap_mat)
y2d <- apply(y3d, c(1, 2), sum)

scr_stan_data <- list(
  n_nonzero_histories = nrow(y2d),
  n_trap = nrow(grid_objs$traplocs), 
  n_occasion = K.jaguar, 
  n_grid = nrow(grid_objs$grid_sf), 
  grid_pts = st_coordinates(grid_objs$grid_sf),
  X = grid_objs$traplocs, 
  y = y2d, 
  n0_prior_scale = 10
)

save(scr_stan_data, file = "./clean_data/scr_stan_data.rda")



