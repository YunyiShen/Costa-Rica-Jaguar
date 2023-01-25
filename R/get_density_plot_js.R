library(rstan)
library(sf)
library(sp)
library(ggplot2)
library(reshape)
source("./R/utils.R")

load("./res/js_stan_fit.rda")
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
  theme_classic() + 
  xlab("Year") +
  ylab("Population size") + 
  geom_point(data = NN_mean) + 
  geom_line(aes(group = 1),data = NN_mean)

ggplot2::ggsave("./res/Figs/js_pop_est.png", width = 6, height = 4, unit = "in")


s <- rstan::extract(m_fit, c("s"))$s

density_est <- list()
png("./res/Figs/js_den_est.png", width = 6 * 2, height = 4 * 2, units = "in",res = 500)
par(mfrow = c(2,2))
for(i in 1:4){
  density_est[[i]] <- JSdensity(s,z,JS_stan_data$grid_pts,i,
                    nx = 40, ny = 40, main = i+2017, 
                    Xl = 0.999*198500/grid_objs$scaling, Xu = 1.001*max(JS_stan_data$grid_pts[,1]),
                    Yl = 0.999*min(JS_stan_data$grid_pts[,2]), 
                    Yu = 1.001*max(JS_stan_data$grid_pts[,2]),
                    )
  points(grid_objs$traplocs, pch = 2)
  points(JS_stan_data$grid_pts, pch = 20, 
         col = adjustcolor("red", alpha.f = 0.2), cex = 0.3)
  plot(as.data.frame(park_boundry$geometry)/grid_objs$scaling, add = T)
}
dev.off()

plot(m_fit, pars = c("beta_env"),plotfun = "trace")
plot(m_fit, pars = c("psi","gamma","phi","p0","alpha1"),plotfun = "trace")
