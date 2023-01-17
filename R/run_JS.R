library(rstan)
options(mc.cores = parallel::detectCores()/2)

load("./clean_data/js_stan_data.rda")
m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp.stan")
m_fit <- sampling(m_init,  data = JS_stan_data,chains = 1, iter = 100, 
        init = function() list(alpha1 = .1, p0 = .01, beta = rep(0,JS_stan_data$n_env)))
#traceplot(m_fit, pars = c("alpha1", "alpha0", "N"))
save(m_fit, file="./res/js_stan_fit.rda")
