library(rstan)
options(mc.cores = parallel::detectCores()/2)

m_init <- stan_model("./stan/scr0-poisson-integrated.stan")

