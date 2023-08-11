## run SCR, with stan
library(rstan)
library(jsonlite)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)
config <- fromJSON("./config.json")

load("./clean_data/Guanecaste/scr_stan_data.rda")

m_init <- stan_model("./stan/discrete-ipp-scr.stan")
scr_stan_data$envX <- 0*scr_stan_data$envX # this line is temporary since the envX is now a place holder
set.seed(12345)
m_fit <- sampling(m_init,  data = scr_stan_data, iter = config$stan_iters, 
                  init = function() list(psi = 0.7, 
                                         alpha1 = 1, 
                                         p0 = .1, beta = rep(0,scr_stan_data$n_env)), 
                  verbose = TRUE)
#m_fit <- vb(m_init,  data = scr_stan_data, iter = config$stan_iters, 
#                  init = function() list(psi = 0.15, 
#                                         alpha1 = .7, 
#                                         p0 = .03, beta = rep(0,scr_stan_data$n_env)))

load("./clean_data/Guanecaste/grid_objs_data.rda")
load("./clean_data/Guanecaste/scr_stan_data.rda")
load("./clean_data/Guanecaste/jaguar_trap_mats_scr.rda")
save(m_fit, grid_objs, jaguar_trap_mats, scr_stan_data, config,
    file="./res/Guanecaste_scr_stan_fit.rda")
