library(rstan)
library(sf)
library(sp)
library(ggplot2)
library(reshape)
library(fileds)
source("./R/utils.R")

load("./res/js_lock_null_stan_fit.rda")
load("./clean_data/grid_objs_data.rda")
load("./clean_data/js_stan_data.rda")
park_boundry <- st_read("./data/mask/Corcovado/Corcovado.shp") |> 
  st_transform(crs = CRS("+proj=utm +zone=17"))

years <- 1:4+2017

z <- rstan::extract(m_fit, c("z"))$z
NN <- apply(z,c(1,3),function(w){sum(w==2)}) |> as.data.frame() # total population size
colnames(NN) <- years
NN <- melt(NN)
NN_mean <- aggregate(value~variable, data = NN, FUN = median)

ggplot(NN, aes(x=variable, y=value)) + 
  geom_violin() + 
  #geom_boxplot() +
  theme_classic() + 
  xlab("Year") +
  ylab("Population size") + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean)

ggplot2::ggsave("./res/Figs/js_pop_est.png", width = 6, height = 4, unit = "in")


s <- rstan::extract(m_fit, c("s"))$s

density_est <- list()
png("./res/Figs/js_null_den_est.png", width = 6 * 2, height = 4 * 2, units = "in",res = 500)
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
            zlim = c(0.1,20), xlab = "", ylab = "", 
            main = i+2017,
            col = gray.colors(20, start = 0., 
                    end = 0.9, gamma = 1.1, rev = TRUE))
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[i,,])>0,], pch = 2)
  for(j in 1:13){ # 13 seen individuals
    points(grid_objs$traplocs[rowSums(JS_stan_data$y [j,,,i])>0,], pch = j+2, col = "blue")
  }
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
}
dev.off()


plot(m_fit, pars = c("beta_env"),plotfun = "hist")
plot(m_fit, pars = c("psi","gamma","phi","p0","alpha1"),plotfun = "trace")
pairs(m_fit, pars = c("beta_env"))

png("./res/Figs/js_act_center.png", width = 10 * 2, height = 4 * 2, units = "in",res = 500)
par(mfrow = c(3,6))
year <- 3
for(i in 1:18){
  JSdensity(s,z,JS_stan_data$grid_pts,year,TRUE,
            nx = 46, ny = 37, main = paste(year+2017,"inid:", i), 
            Xl = min(JS_stan_data$grid_pts[,1])-.2, 
            Xu = max(JS_stan_data$grid_pts[,1])+.2,
            Yl = min(JS_stan_data$grid_pts[,2])-.2, 
            Yu = max(JS_stan_data$grid_pts[,2])+.2,whichguy = i
  )
  points(grid_objs$traplocs[rowSums(JS_stan_data$deploy[year,,])>0,], pch = 2)
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  points(JS_stan_data$grid_pts[s[,i],],col = "blue")
  for(j in 1:4){
    points(grid_objs$traplocs[rowSums(JS_stan_data$y[i,,,j])>0,], pch = 10, col = "purple")
  }
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
}
dev.off()
