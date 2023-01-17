library(rstan)
options(mc.cores = parallel::detectCores()/2)

load("./clean_data/js_stan_data.rda")
m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp-binomial.stan")
set.seed(42)
m_fit <- sampling(m_init,  data = JS_stan_data,chains = 1, iter = 1000, 
        init = function() list(alpha1 = .1, p0 = .2, beta = rep(0,JS_stan_data$n_env)))
#traceplot(m_fit, pars = c("alpha1", "alpha0", "N"))
save(m_fit, file="./res/js_stan_fit.rda")
