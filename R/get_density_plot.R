library(rstan)
library(sf)
source("./R/utils.R")

load("./res/scr_stan_fit.rda")
load("./clean_data/grid_objs_data.rda")

stan_plot(m_fit, pars = c("alpha1", "alpha0", "N"))
stan_dens(m_fit, pars = c("alpha1", "alpha0", "N"))

SCR0density(m_fit)
points(grid_objs$traplocs, pch = 2)
for(i in 1:15)
  points(grid_objs$traplocs[jaguar_trap_mats$cap_mat$trap[jaguar_trap_mats$cap_mat$ind_id==i],], col = "blue", pch = i)

plot_s(3, m_fit, add = T)
plot_s(6, m_fit, add = T)
plot_s(8, m_fit, add = T)
plot_s(9, m_fit, add = T)
