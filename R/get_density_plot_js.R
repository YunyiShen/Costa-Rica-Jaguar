library(rstan)
library(sf)
library(sp)
library(ggplot2)
library(reshape)
library(fields)
source("./R/utils.R")
source("./R/forecast_js.R")

load("./res/js_lock_stan_fit.rda")
#load("./clean_data/grid_objs_data.rda")
#load("./clean_data/js_stan_data.rda")
park_boundry <- st_read("./data/mask/Corcovado/Corcovado.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

old_trap_polygon <- st_read("./data/2007_map/trap_polygon_buffer.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))
old_trap_polygon_inner <- st_read("./data/2007_map/trap_polygon.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

years <- 1:7+2014
year_out <- 2
years_w_forcast <- 1:(7+year_out) + 2014 

#### pop size ####
area_size <- nrow(JS_stan_data$grid_pts) # we have 1km^2 grids

z <- rstan::extract(m_fit, c("z"))$z
forecast_z <- forecast_js(m_fit, year_out, recruit_factor = 1)
NN <- apply(z,c(1,3),function(w){sum(w==2)}) |> as.data.frame() # total population size
NN <- cbind(NN, apply(forecast_z,c(1,3),function(w){sum(w==2)}) |> as.data.frame())
colnames(NN) <- years_w_forcast
NN <- melt(NN)
NN_mean <- aggregate(value~variable, data = NN, FUN = median)

jpeg("./res/Figs/js_prey_avg_den_est.jpg", width = 6, height = 4, units = "in",res = 500)
ggplot(NN, aes(x=variable, y=value/area_size * 100)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("") +
  #ylab(paste0("Density (/100km",expression(km^2),")")) + 
  ylab(expression(Density ~ group("(",100~km^2,")")))+
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean) + 
  geom_vline(xintercept = 7.5, lty = 2) + 
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line())
  #geom_rect(aes(xmin=7.5, xmax=9.5, ymin=0, ymax=3), alpha = .0002) 
dev.off()

#ggplot2::ggsave("./res/Figs/js_prey_avg_den_est.jpg", width = 6, height = 4, unit = "in")

jpeg("./res/Figs/js_prey_pop_est.jpg", width = 6, height = 4, units = "in",res = 500)
ggplot(NN, aes(x=variable, y=value)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("") +
  ylab("Population size") + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean)+ 
  geom_vline(xintercept = 7.5, lty = 2) + 
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line())
dev.off()
#ggplot2::ggsave("./res/Figs/js_prey_pop_est.jpg", width = 6, height = 4, unit = "in")


#### local density ####
s <- rstan::extract(m_fit, c("s"))$s
density_est <- list()
jpeg("./res/Figs/js_prey_den_est.jpg", width = 7.5 * 1.2, height = 4 * 1.2, units = "in",res = 500)
par(mfrow = c(2,4),mar = c(.5,.5,1.7,.5), mgp = c(1.5, 0.5, 0))
for(i in 1:7){
  density_est[[i]] <- JSdensity(s,z,JS_stan_data$grid_pts,i,TRUE,
                    nx = 46, ny = 37, main = i+2014, 
                    Xl = min(JS_stan_data$grid_pts[,1])-.2, 
                    Xu = max(JS_stan_data$grid_pts[,1])+.2,
                    Yl = min(JS_stan_data$grid_pts[,2])-.2, 
                    Yu = max(JS_stan_data$grid_pts[,2])+.2,
                    plotit = FALSE
                    )
  image(density_est[[i]]$grid$xg, 
        density_est[[i]]$grid$yg, 
        density_est[[i]]$Dn, 
        zlim = c(0.,30), xlab = "", ylab = "", 
        main = i+2014,
        col = gray.colors(30, start = 0., 
                          end = 0.9, gamma = .4, rev = TRUE), xaxt='n', yaxt='n',font.main = 1)
  
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[,,i])>0,], pch = 1)
  for(j in 1:13){ # 13 seen individuals
    points(grid_objs$traplocs[rowSums(JS_stan_data$y [j,,,i])>0,], pch = 19)
  }
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  plot(as.data.frame(old_trap_polygon$geometry)/grid_objs$scaling, add = T, lty = 2)
  plot(as.data.frame(old_trap_polygon_inner$geometry)/grid_objs$scaling, add = T, lty = 2)
  if(i==1){
    #legend("topright", legend = c("Camera traps",
    #                              "Jaguars detected",
    #                              "Grid points",
    #                              "Salom-Pérez et al. (2007)"),
    #       pch = c(1,19,20,NA), cex = c(1,1,1,1), 
    #       lty = c(NA,NA,NA,2),
    #       col = c("black","black",adjustcolor("red", alpha.f = 0.2),"black"))
  }
}

#fields::image.plot(density_est[[i]]$grid$xg, 
#                   density_est[[i]]$grid$yg, 
#                   density_est[[i]]$Dn, 
#                   zlim = c(0.,30), xlab = "", ylab = "", 
#                   main = i+2014,
#                   col = gray.colors(30, start = 0., 
#                                     end = 0.9, gamma = .4, rev = TRUE), 
#                   legend.mar = -60, legend.only = T)

image(density_est[[i]]$grid$xg, 
      density_est[[i]]$grid$yg, 
      density_est[[i]]$Dn + NA, 
      zlim = c(0.,30), xlab = "", ylab = "", 
      main = "",
      col = gray.colors(30, start = 0., 
                        end = 0.9, gamma = .4, rev = TRUE), 
      bty ="n",
      xaxt='n', yaxt='n',font.main = 1)

legend("topright", 
       #inset=c(1.2, 0),
       legend = c("Camera traps",
                              "Jaguars detected",
                              "Grid points",
                              "Salom-Pérez et al.\n(2007)",
                              "Corcovado\nNational Park"
                  ),
       pch = c(1,19,20,NA,NA), cex = c(1,1,1,1,1), 
       lty = c(NA,NA,NA,2,NA),
       fill=c(NA,NA,NA,NA,"gray80"),
       border = c(NA,NA,NA,NA,"black"),
       bty = "n",
       col = c("black","black",adjustcolor("red", alpha.f = 0.2),
               "black"),
       y.intersp = 2, 
       xjust = 0
       )

fields::image.plot(density_est[[i]]$grid$xg, 
                   density_est[[i]]$grid$yg, 
                   density_est[[i]]$Dn + NA, 
                   zlim = c(0.,30), xlab = "", ylab = "", 
                   main = "",
                   col = gray.colors(30, start = 0., 
                                     end = 0.9, gamma = .4, rev = TRUE), 
                   axes = FALSE,    # Suppress axes
                   box = FALSE,     # Remove bounding box
                   legend.mar = 54.535, 
                   legend.line = 1,
                   #legend.lab = "Density",
                   legend.only = T, 
                   legend.args = list(text = "Density",cex = 0.7,adj = -0.005, line = 0.3)
                   )

#legend("center", 
#       inset=c(1.2, 0),
#       legend = c("Camera traps",
#                              "Jaguars detected",
#                              "Grid points",
#                              "Salom-Pérez et al. (2007)"),
#       pch = c(1,19,20,NA), cex = c(1,1,1,1), 
#       lty = c(NA,NA,NA,2),
#       col = c("black","black",adjustcolor("red", alpha.f = 0.2),"black", xpd=TRUE))

dev.off()




#### sd ####
jpeg("./res/Figs/js_prey_den_est_relative_sd.jpg", width = 8.5 * 1.5, height = 4 * 1.5, units = "in",res = 500)
par(mfrow = c(2,4),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
for(i in 1:7){
  density_est[[i]] <- JSdensity(s,z,JS_stan_data$grid_pts,i,TRUE,
                                nx = 46, ny = 37, main = i+2014, 
                                Xl = min(JS_stan_data$grid_pts[,1])-.2, 
                                Xu = max(JS_stan_data$grid_pts[,1])+.2,
                                Yl = min(JS_stan_data$grid_pts[,2])-.2, 
                                Yu = max(JS_stan_data$grid_pts[,2])+.2,
                                plotit = FALSE, calculate_sd = TRUE
  )
  fields::image.plot(density_est[[i]]$grid$xg, 
                     density_est[[i]]$grid$yg, 
                     (density_est[[i]]$Dn_sd+1/2000)/(density_est[[i]]$Dn+1/2000), 
                     zlim = c(0.,50), xlab = "", ylab = "", 
                     main = i+2014,
                     col = gray.colors(30, start = 0., 
                                       end = 0.9, gamma = .8, rev = TRUE), legend.mar = 7)
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[,,i])>0,], pch = 2)
  for(j in 1:13){ # 13 seen individuals
    points(grid_objs$traplocs[rowSums(JS_stan_data$y [j,,,i])>0,], pch = j+2, col = "blue")
  }
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  plot(as.data.frame(old_trap_polygon$geometry)/grid_objs$scaling, add = T, lty = 2)
  plot(as.data.frame(old_trap_polygon_inner$geometry)/grid_objs$scaling, add = T, lty = 2)
  if(i==1){
    legend("topright", legend = c("deployed traps",
                                  "seen individuals",
                                  "grid points in study area",
                                  "Salom-Perez et al. 2007"),
           pch = c(2,3,20,NA), cex = c(1,1,1,1), 
           lty = c(NA,NA,NA,2),
           col = c("black","blue",adjustcolor("red", alpha.f = 0.2),"black"))
  }
}
dev.off()



#### some sanity checks ####
#load("./clean_data/jaguar_trap_mats_js.rda")
jaguar_id <- jaguar_trap_mats$ids$ind_ids
inid_plot <- 1:nrow(jaguar_id)

jpeg("./res/Figs/js_prey_act_center.jpg", width = 10 * 2, height = 6 * 2, units = "in",res = 500)


par(mfrow = c(3,4),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
year <- 7
avg_in <- 0 * inid_plot

for(i in inid_plot){
  si <- s[,i]
  s_loc <- JS_stan_data$grid_pts[si, ] * grid_objs$scaling 
  s_loc <- lapply(1:nrow(s_loc), function(i,s_loc){
    st_point(s_loc[i,])
  }, s_loc) |>
    st_sfc(crs = "+proj=utm +zone=17")
  lst <- st_intersects( park_boundry,s_loc)[[1]]  |> length()
  avg_in[i] <- lst/nrow(s)
  
  JSdensity(s,z,JS_stan_data$grid_pts,year,TRUE,
            nx = 46, ny = 37, main = jaguar_id$jaguar[i], 
            Xl = min(JS_stan_data$grid_pts[,1])-.2, 
            Xu = max(JS_stan_data$grid_pts[,1])+.2,
            Yl = min(JS_stan_data$grid_pts[,2])-.2, 
            Yu = max(JS_stan_data$grid_pts[,2])+.2, xlab = "", ylab = "", plot_scale = (i%%4 == 0)
  )
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[year,,])>0,], pch = 2)
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  points(JS_stan_data$grid_pts[s[,i],],col = adjustcolor("blue", alpha.f = 0.2), pch = 1, cex = 0.5)
  
    
  for(j in 1:7){
    points(grid_objs$traplocs[rowSums(JS_stan_data$y[i,,,j])>0,], pch = 17, col = "purple", cex = 1.2)
  }
  
  sort_id <- sort(s[,i])
  #mode_idx <- as.numeric(names(which.max(summary( as.factor( s[,i] )))))
  mode_idx <- sort_id[which(!duplicated(sort(s[,i])))[which.max( diff(which(!duplicated(sort(s[,i]))) )) ]]
  
  points(JS_stan_data$grid_pts[mode_idx,1],  JS_stan_data$grid_pts[mode_idx,2],
         col = "gold", pch = 16, cex = 1.2)
  
  if(i==1){
    legend("topright", legend = c("deployed traps",
                                  "trap with detections",
                                  "posterior of activities center",
                                  "posterior mode of activities center",
                                  "grid points in study area"),
           pch = c(2,17,1,16,20), cex = c(1,1,1,1,1),col = c("black","purple","blue","gold",adjustcolor("red", alpha.f = 0.2)))
  }
  text(20.55, 93.5, label = paste0(signif(avg_in[i],3)*100,"% in"))
}
dev.off()

#### parameters####
params <- rstan::extract(m_fit, c("psi","gamma","phi","p0","alpha1","beta_env"))
jpeg("./res/Figs/js_prey_vital_rate.jpg", width = 6 * 1.2, height = 4 * 1.2, units = "in",res = 500)
par(mfrow = c(2,3),mar = c(2.5,2.5,1.7,.5), mgp = c(1.5, 0.5, 0))

plot(density(params$psi*60), main = "Initial recruitment", xlab = "",font.main = 1)
polygon(density(params$psi*60), col = "#9b9b9b")

plot(density(params$gamma*60), main = "Annual recruitment", xlab = "",font.main = 1)
polygon(density(params$gamma*60), col = "#9b9b9b")

plot(density(params$phi), main = "Annual survival", xlab = "",font.main = 1)
polygon(density(params$phi), col = "#9b9b9b")

plot(density(params$p0), main = "Baseline detection rate", xlab = "",font.main = 1)
polygon(density(params$p0), col = "#9b9b9b")

plot(density(params$alpha1), main = "Detection rate decay", xlab = "",font.main = 1)
polygon(density(params$alpha1), col = "#9b9b9b")

plot(density(params$beta_env), main = "Effect of prey", xlab = "", xlim = c(-5,2),font.main = 1)
polygon(density(params$beta_env), col = "#9b9b9b")
abline(v=0, lwd = 2)
dev.off()


lw_ci <- sapply(params, quantile, .025)
hi_ci <- sapply(params, quantile, .975)
mi_ci <- sapply(params, quantile, .5)

write.csv(data.frame(median = mi_ci,
  ci_low = lw_ci, 
  ci_high = hi_ci),"./res/js_lock_prey_ci_vital.csv")  


plot(m_fit, pars = c("psi", "gamma", "phi", "p0", "alpha1", "beta_env"), plotfun = "trace")
ggsave("./res/Figs/js_prey_trace.jpg", width = 6, height = 4, unit = "in")

#### density in old area ####
density_sample <- list()
for(i in 1:7){
  density_sample[[years[i]]] <- density_within_polygon(JS_stan_data$grid_pts * grid_objs$scaling,
                                                     old_trap_polygon, s, z, i, locked = TRUE)
}
density_sample <- Reduce(cbind, density_sample) |> 
  as.data.frame()
colnames(density_sample) <- years

density_sample <- melt(density_sample)
density_sample_mean <- aggregate(value~variable, data = density_sample, FUN = median)

ggplot(density_sample, aes(x=variable, y=value)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("Year") +
  ylab("Density (/100km^2) \n in Salom-Pérez et al. 2007 range") + 
  geom_point(data = density_sample_mean) + 
  geom_line(aes(group = 1),data = density_sample_mean)

ggplot2::ggsave("./res/Figs/js_prey_den_est_SP07.jpg", width = 6, height = 4, unit = "in")


#### density in park ####

density_sample <- list()
for(i in 1:7){
  density_sample[[years[i]]] <- density_within_polygon(JS_stan_data$grid_pts * grid_objs$scaling,
                                                       park_boundry, s, z, i, locked = TRUE)
}
density_sample <- Reduce(cbind, density_sample) |> 
  as.data.frame()
colnames(density_sample) <- years

density_sample <- melt(density_sample)
density_sample_mean <- aggregate(value~variable, data = density_sample, FUN = median)

js_prey_den_est_park <- ggplot(density_sample, aes(x=variable, y=value)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("Year") +
  ylab("Density (/100km^2) \n within park boundary") + 
  geom_point(data = density_sample_mean) + 
  geom_line(aes(group = 1),data = density_sample_mean)

js_prey_den_est_park
save(js_prey_den_est_park,file = "./res/Figs/js_prey_den_est_park.rda")

ggplot2::ggsave( "./res/Figs/js_prey_den_est_park.jpg", width = 6, height = 4, unit = "in")

