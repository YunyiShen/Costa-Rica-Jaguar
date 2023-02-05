library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

load("./clean_data/js_stan_data.rda")
m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp-binomial-lock.stan")

envX_static <- JS_stan_data$envX[1,,-1] 
envX_dyna <- JS_stan_data$envX[,,1]|> 
  apply(1, function(w){(w-mean(w))/sd(w)})
envX_dyna <- rowMeans(envX_dyna)
JS_stan_data$envX <- cbind(envX_dyna, envX_static)

set.seed(42)
m_fit <- sampling(m_init,  data = JS_stan_data,chains = 1, iter = 10, 
                  init = function() list(gamma = 0.2, psi = 0.15, 
                                         phi = 0.85, alpha1 = .7, 
                                         p0 = .03, beta = rep(0,JS_stan_data$n_env)), 
                  verbose = TRUE)

