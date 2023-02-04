library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

load("./clean_data/js_stan_data.rda")
m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp-binomial-lock.stan")
