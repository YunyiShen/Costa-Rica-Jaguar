library(rstan)
library(sf)
library(sp)
library(reshape)
options(mc.cores = 4)
rstan_options(auto_write = TRUE)
source("./R/utils.R")
source("./R/forecast_js.R")
library(jsonlite)
config <- fromJSON("./config.json")


#load("./res/js_lock_stan_fit.rda")
#rm(m_fit)
load("./clean_data/grid_objs_data.rda")
load("./clean_data/js_stan_data.rda")
load("./clean_data/jaguar_trap_mats_js.rda")
#### set up the data with only the old range
old_trap_polygon <- st_read("./data/2007_map/trap_polygon_buffer.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

old_trap_polygon_inner <- st_read("./data/2007_map/trap_polygon.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

pts.sf <- st_as_sf(data.frame(x = JS_stan_data$grid_pts[,1] * grid_objs$scaling, 
                                y = JS_stan_data$grid_pts[,2] * grid_objs$scaling, 
                                id = 1:nrow(JS_stan_data$grid_pts)), 
                     coords = c("x","y"), crs =  CRS("+proj=utm +zone=17")) 
good_pts <- st_filter(pts.sf, old_trap_polygon)$id


traps.sf <- st_as_sf(data.frame(x = JS_stan_data$X[,1] * grid_objs$scaling, 
                              y = JS_stan_data$X[,2] * grid_objs$scaling, 
                              id = 1:nrow(JS_stan_data$X)), 
                   coords = c("x","y"), crs =  CRS("+proj=utm +zone=17")) 
good_traps <- st_filter(traps.sf, old_trap_polygon)$id


envX_dyna <- as.matrix(JS_stan_data$envX[,,1])|> 
  apply(2, function(w){(w-mean(w))/sd(w)})
envX_dyna <- rowMeans(envX_dyna)
JS_stan_data$envX <- as.matrix(envX_dyna)


JS_stan_data$n_trap <- length(good_traps)
JS_stan_data$y <- JS_stan_data$y[,good_traps,,]
JS_stan_data$deploy <- JS_stan_data$deploy[good_traps,,]
JS_stan_data$X <- JS_stan_data$X[good_traps,]

JS_stan_data$n_grid <- length(good_pts)
JS_stan_data$grid_pts <- JS_stan_data$grid_pts[good_pts,]

JS_stan_data$envX <- as.matrix(JS_stan_data$envX[good_pts,])

m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp-binomial-lock.stan")

set.seed(42)
m_fit <- sampling(m_init,  data = JS_stan_data,chains = 4, iter = 1000, 
                  init = function() list(gamma = 0.2, psi = 0.15, 
                                         phi = 0.85, alpha1 = .7, 
                                         p0 = .03, beta = rep(0,JS_stan_data$n_env)), 
                  verbose = TRUE)
save(m_fit, grid_objs, jaguar_trap_mats, JS_stan_data, config, file="./res/js_lock_stan_fit_07_range.rda")

#m_fit <- vb(m_init,  data = JS_stan_data, 
#                  init = function() list(gamma = 0.2, psi = 0.15, 
#                                         phi = 0.85, alpha1 = .7, 
#                                         p0 = .03, beta = rep(0,JS_stan_data$n_env)))

#### get average density ####
load("./res/js_lock_stan_fit_07_range.rda")
area_size <- nrow(JS_stan_data$grid_pts) # we have 1km^2 grids
z <- rstan::extract(m_fit, c("z"))$z
years <- 1:7+2014
year_out <- 0
years_w_forcast <- 1:(7+year_out) + 2014 

forecast_z <- forecast_js(m_fit, year_out, recruit_factor = 1)
NN <- apply(z,c(1,3),function(w){sum(w==2)}) |> as.data.frame() # total population size
#NN <- cbind(NN, apply(forecast_z,c(1,3),function(w){sum(w==2)}) |> as.data.frame())
colnames(NN) <- years_w_forcast
NN <- melt(NN)
NN_mean <- aggregate(value~variable, data = NN, FUN = median)

png("./res/Figs/js_prey_avg_den_est_restricted_area.png", width = 6, height = 4, units = "in",res = 500)
js_prey_avg_den_est_restricted_area <- ggplot(data = NN, aes(x=variable, y=value/area_size * 100))  + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("") +
  ylab(expression(Density ~ group("(",100~km^2,")")))+
  geom_rect(aes(xmin=0.45, xmax=7.7+year_out, ymin=6.98-2.36, ymax=6.98+2.36), alpha = .002) + 
  geom_violin() + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean) + 
  geom_hline(yintercept = 6.98, show.legend = "Salom-Pérez et al. 2007",
             linetype=2) + 
  geom_label(x = 0.99, y = 9.25, label = "Salom-Pérez\n et al. 2007",
             #color="", 
             size=2.5 ) + 
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line()) #+ 
  #geom_vline(xintercept = 7.5, lty = 2)
dev.off()
#ggplot2::ggsave("./res/Figs/js_prey_avg_den_est_restricted_area.png", width = 6, height = 4, unit = "in")


### combine with Fig2 
load("./res/Figs/js_prey_den_est_park.rda")
js_prey_den_est_park <- js_prey_den_est_park + 
  theme(panel.grid.major.y = element_line(),
        panel.grid.minor.y = element_line(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  xlab("") + 
  ylab(expression(Density ~ "in" ~ park ~ group("(",100~km^2,")"))) 
  
library(ggpubr)
jpeg("./res/Figs/js_prey_avg_den_est_restricted_area_comb.jpg", width = 5.5/1.15, height = 4.5/1., unit = "in", res = 500)
ggarrange(js_prey_den_est_park + theme(plot.margin = unit(c(0,0.1,0,0), 'lines')), 
          js_prey_avg_den_est_restricted_area + theme(plot.margin = unit(c(0,0.1,0,0), 'lines')), 
          nrow = 2, labels=c("(a)","(b)"), align = "hv", hjust = -3.,font.label = list(face = "plain"))
dev.off()


#### plot density ####
z <- rstan::extract(m_fit, c("z"))$z
s <- rstan::extract(m_fit, c("s"))$s
density_est <- list()
png("./res/Figs/js_prey_den_est_restricted_map.png", width = 8.5 * 1.5, height = 4 * 1.5, units = "in",res = 500)
par(mfrow = c(2,4),mar = c(2.5,2.5,1,.5), mgp = c(1.5, 0.5, 0))
for(i in 1:7){
  density_est[[i]] <- JSdensity(s,z,JS_stan_data$grid_pts,i,TRUE,
                                nx = 15, ny = 18, main = i+2014, 
                                Xl = min(JS_stan_data$grid_pts[,1])-.2, 
                                Xu = max(JS_stan_data$grid_pts[,1])+.2,
                                Yl = min(JS_stan_data$grid_pts[,2])-.2, 
                                Yu = max(JS_stan_data$grid_pts[,2])+.2,
                                plotit = FALSE
  )
  fields::image.plot(density_est[[i]]$grid$xg, 
                     density_est[[i]]$grid$yg, 
                     density_est[[i]]$Dn, 
                     zlim = c(0.,140), xlab = "", ylab = "", 
                     main = i+2014,
                     col = gray.colors(140, start = 0., 
                                       end = 0.9, gamma = .6, rev = TRUE), 
                     legend.mar = 8.5)
  points(grid_objs$traplocs[good_traps,][rowSums(JS_stan_data$deploy[i,,])>0,], pch = 2)
  for(j in 1:13){ # 13 seen individuals
    points(grid_objs$traplocs[good_traps,][rowSums(JS_stan_data$y [j,,,i])>0,], pch = j+2, col = "blue")
  }
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  #plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
  plot(as.data.frame(old_trap_polygon$geometry)/grid_objs$scaling, add = T, lty = 2)
  plot(as.data.frame(old_trap_polygon_inner$geometry)/grid_objs$scaling, add = T, lty = 2)
  if(i==10){
    legend("topright", legend = c("deployed traps",
                                  "seen individuals",
                                  "grid points in study area",
                                  "Salom-Pérez et al. 2007"),
           pch = c(2,3,20,NA), cex = c(1,1,1,1), 
           lty = c(NA,NA,NA,2),
           col = c("black","blue",adjustcolor("red", alpha.f = 0.2),"black"))
  }
}
dev.off()




