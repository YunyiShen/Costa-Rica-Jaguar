library(rstan)
options(mc.cores = parallel::detectCores()/2)
rstan_options(auto_write = TRUE)

load("./clean_data/js_stan_data.rda")
m_init <- stan_model("./stan/spatial-jolly-seber-discrete-ipp-binomial.stan")
set.seed(42)
m_fit <- sampling(m_init,  data = JS_stan_data,chains = 2, iter = 1000, 
        init = function() list(gamma = 0.2, psi = 0.15, 
                                phi = 0.85, alpha1 = .7, 
                                p0 = .03, beta = rep(0,JS_stan_data$n_env)), 
        verbose = TRUE)

#m_fit <- rstan::optimizing(m_init,  data = JS_stan_data)
#traceplot(m_fit, pars = c("alpha1", "alpha0", "N"))
save(m_fit, file="./res/js_stan_fit.rda")

# z <- matrix(1, JS_stan_data$M, JS_stan_data$T)
# for(i in 1:JS_stan_data$M){
#     for(tt in 1:JS_stan_data$T){
#         z[i,tt] <- `[`(m_fit$par, paste0("z[",i,",",tt,"]"))
#     } 
# }

# s <- matrix(1, JS_stan_data$M, JS_stan_data$T)
# for(i in 1:JS_stan_data$M){
#     for(tt in 1:JS_stan_data$T){
#         s[i,tt] <- `[`(m_fit$par, paste0("s[",i,",",tt,"]"))
#     } 
# }

# plot(JS_stan_data$grid_pts)
# for(i in 1:JS_stan_data$M){
#    for(tt in 1:JS_stan_data$T){
#        if(z[i,tt] == 2){
#            points(JS_stan_data$grid_pts[s[i,tt],1],
#            JS_stan_data$grid_pts[s[i,tt],2], col = "red")
#        }
#     }
# }

