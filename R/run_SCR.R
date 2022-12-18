## run SCR, with stan
library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

load("./clean_data/scr_stan_data.rda")

m_init <- stan_model("./stan/scr0-poisson-integrated.stan")
m_fit <- sampling(m_init,  data = scr_stan_data)
traceplot(m_fit, pars = c("alpha1", "alpha0", "N"))

save(m_fit, file="./res/scr_stan_fit.rda")