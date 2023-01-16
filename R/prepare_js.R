library(sp)
library(sf)
library(terra)
source("./R/utils.R")

load("./clean_data/jaguar_trap_mats_js.rda")
load("./clean_data/grid_objs_data.rda")

M <- 50 # the augmentation population size



