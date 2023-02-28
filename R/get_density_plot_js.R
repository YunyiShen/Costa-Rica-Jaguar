library(rstan)
library(sf)
library(sp)
library(ggplot2)
library(reshape)
library(fields)
source("./R/utils.R")

load("./res/js_lock_prey_stan_fit.rda")
load("./clean_data/grid_objs_data.rda")
load("./clean_data/js_stan_data.rda")
park_boundry <- st_read("./data/mask/Corcovado/Corcovado.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

old_trap_polygon <- st_read("./data/2007_map/trap_polygon_buffer.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))
old_trap_polygon_inner <- st_read("./data/2007_map/trap_polygon.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

years <- 1:4+2017

#### pop size ####
area_size <- nrow(JS_stan_data$grid_pts) # we have 1km^2 grids

z <- rstan::extract(m_fit, c("z"))$z
NN <- apply(z,c(1,3),function(w){sum(w==2)}) |> as.data.frame() # total population size
colnames(NN) <- years
NN <- melt(NN)
NN_mean <- aggregate(value~variable, data = NN, FUN = median)

ggplot(NN, aes(x=variable, y=value/area_size * 100)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("Year") +
  ylab("Density (/100km^2)") + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean)

ggplot2::ggsave("./res/Figs/js_prey_avg_den_est.png", width = 6, height = 4, unit = "in")

ggplot(NN, aes(x=variable, y=value)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("Year") +
  ylab("Population size") + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean)

ggplot2::ggsave("./res/Figs/js_prey_pop_est.png", width = 6, height = 4, unit = "in")


#### local density ####
s <- rstan::extract(m_fit, c("s"))$s
density_est <- list()
png("./res/Figs/js_prey_den_est.png", width = 6 * 1.5, height = 4 * 1.5, units = "in",res = 500)
par(mfrow = c(2,2),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
for(i in 1:4){
  density_est[[i]] <- JSdensity(s,z,JS_stan_data$grid_pts,i,TRUE,
                    nx = 46, ny = 37, main = i+2017, 
                    Xl = min(JS_stan_data$grid_pts[,1])-.2, 
                    Xu = max(JS_stan_data$grid_pts[,1])+.2,
                    Yl = min(JS_stan_data$grid_pts[,2])-.2, 
                    Yu = max(JS_stan_data$grid_pts[,2])+.2,
                    plotit = FALSE
                    )
  fields::image.plot(density_est[[i]]$grid$xg, 
            density_est[[i]]$grid$yg, 
            density_est[[i]]$Dn, 
            zlim = c(0.,18), xlab = "", ylab = "", 
            main = i+2017,
            col = gray.colors(20, start = 0., 
                    end = 0.9, gamma = .8, rev = TRUE))
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[i,,])>0,], pch = 2)
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
load("./clean_data/jaguar_trap_mats_js.rda")
jaguar_id <- jaguar_trap_mats$ids$ind_ids
inid_plot <- 1:nrow(jaguar_id)

png("./res/Figs/js_prey_act_center.png", width = 10 * 2, height = 6 * 2, units = "in",res = 500)
par(mfrow = c(3,4),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
year <- 4
for(i in inid_plot){
  JSdensity(s,z,JS_stan_data$grid_pts,year,TRUE,
            nx = 46, ny = 37, main = jaguar_id$jaguar[i], 
            Xl = min(JS_stan_data$grid_pts[,1])-.2, 
            Xu = max(JS_stan_data$grid_pts[,1])+.2,
            Yl = min(JS_stan_data$grid_pts[,2])-.2, 
            Yu = max(JS_stan_data$grid_pts[,2])+.2, xlab = "", ylab = ""
  )
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[year,,])>0,], pch = 2)
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  points(JS_stan_data$grid_pts[s[,i],],col = "blue")
  for(j in 1:4){
    points(grid_objs$traplocs[rowSums(JS_stan_data$y[i,,,j])>0,], pch = 10, col = "purple")
  }
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  if(i==1){
    legend("topright", legend = c("deployed traps",
                                  "trap with detections",
                                  "posterior of activities center",
                                  "grid points in study area"),
           pch = c(2,10,1,20), cex = c(1,1,1,1),col = c("black","purple","blue",adjustcolor("red", alpha.f = 0.2)))
  }
  
}
dev.off()

#### parameters####
params <- rstan::extract(m_fit, c("psi","gamma","phi","p0","alpha1","beta_env"))
png("./res/Figs/js_prey_vital_rate.png", width = 6 * 1.5, height = 4 * 1.5, units = "in",res = 500)
par(mfrow = c(2,3),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))

plot(density(params$psi*60), main = "initial recruitment", xlab = "")
polygon(density(params$psi*60), col = "#9b9b9b")

plot(density(params$gamma*60), main = "annual recruitment", xlab = "")
polygon(density(params$gamma*60), col = "#9b9b9b")

plot(density(params$phi), main = "annual survival", xlab = "")
polygon(density(params$phi), col = "#9b9b9b")

plot(density(params$p0), main = "baseline detection rate", xlab = "")
polygon(density(params$p0), col = "#9b9b9b")

plot(density(params$alpha1), main = "detection rate decay", xlab = "")
polygon(density(params$alpha1), col = "#9b9b9b")

plot(density(params$beta_env), main = "effect of prey", xlab = "")
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
ggsave("./res/Figs/js_prey_trace.png", width = 6, height = 4, unit = "in")

#### density in old area ####
density_sample <- list()
for(i in 1:4){
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
  ylab("Density (/100km^2) \n in Salom-Perez et al. 2007 range") + 
  geom_point(data = density_sample_mean) + 
  geom_line(aes(group = 1),data = density_sample_mean)

ggplot2::ggsave("./res/Figs/js_prey_den_est_SP07.png", width = 6, height = 4, unit = "in")


#### density in park ####

density_sample <- list()
for(i in 1:4){
  density_sample[[years[i]]] <- density_within_polygon(JS_stan_data$grid_pts * grid_objs$scaling,
                                                       park_boundry, s, z, i, locked = TRUE)
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
  ylab("Density (/100km^2) \n within park boundary") + 
  geom_point(data = density_sample_mean) + 
  geom_line(aes(group = 1),data = density_sample_mean)

ggplot2::ggsave("./res/Figs/js_prey_den_est_park.png", width = 6, height = 4, unit = "in")

