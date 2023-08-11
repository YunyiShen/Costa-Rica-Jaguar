library(rstan)
library(sf)
source("./R/utils.R")

load("./res/Guanecaste_scr_stan_fit.rda")
#load("./clean_data/Guanecaste/grid_objs_data.rda")
#load("./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")

park_boundry <- st_read("./data/Guanecaste/shore_line/CRI_adm0.shp") |> 
  st_transform(crs = sp::CRS("+proj=utm +zone=17"))

stan_trace(m_fit, pars = c("p0", "alpha1", "psi", "beta_env"))
stan_dens(m_fit, pars = c("p0", "alpha1","psi", "beta_env"))

z <- rstan::extract(m_fit, c("z"))$z
NN <- rowSums((z))
s <- rstan::extract(m_fit, c("s"))$s

png("./res/Figs_Guanecaste/scr_vanilla_est.png", width = 5 , height = 3 , 
  units = "in",res = 500)
par(mar = c(2.5,2.5,1,.1), mgp = c(1.5, 0.5, 0))
density_est<- SCRdensity(s,z,scr_stan_data$grid_pts,TRUE,
                    nx = 44, ny = 33, 
                    Xl = min(scr_stan_data$grid_pts[,1])-.2, 
                    Xu = max(scr_stan_data$grid_pts[,1])+.2,
                    Yl = min(scr_stan_data$grid_pts[,2])-.2, 
                    Yu = max(scr_stan_data$grid_pts[,2])+.2,
                    plotit = FALSE
                    )
  fields::image.plot(density_est$grid$xg, 
            density_est$grid$yg, 
            density_est$Dn, 
            zlim = c(0.,100), xlab = "", ylab = "", 
            col = gray.colors(30, start = 0., 
                    end = 0.9, gamma = .4, rev = TRUE), legend.mar = 7)
  points(grid_objs$traplocs[rowSums(scr_stan_data$deploy[,])>0,], pch = 2)
  points(scr_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  legend("bottomleft", legend = c("deployed traps",
                                  "grid points in study area"),
           pch = c(2,20), cex = c(1,1),
           col = c("black",adjustcolor("red", alpha.f = 0.2)))
dev.off()
