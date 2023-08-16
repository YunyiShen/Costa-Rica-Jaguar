forecast_js <- function(m_fit, year_out = 2, recruit_factor = 1){
  useful_var <- extract(m_fit, c("z","gamma","phi"))
  useful_var$z <- useful_var$z[,,dim(useful_var$z)[3]]
  
  z_new <- array(NA, dim = c(nrow(useful_var$z), ncol(useful_var$z), year_out + 1))
  z_new[,,1] <- useful_var$z
  useful_var$z <- NULL
  
  for(tt in 1:year_out + 1){
    survive_or_not <- matrix(runif(prod(dim(z_new[,,1]))), 
                             nrow = nrow(z_new[,,1])) <= 
      (useful_var$phi %*% t(rep(1, ncol(z_new[,,1]))))
    
    recruit_or_not <- matrix(runif(prod(dim(z_new[,,1]))), 
                             nrow = nrow(z_new[,,1])) <= 
      recruit_factor * (useful_var$gamma %*% t(rep(1, ncol(z_new[,,1]))))
    
    tmp <- z_new[,,tt-1]
    tmp[ tmp==2 & !survive_or_not ] <- 3
    tmp[tmp == 1 & recruit_or_not] <- 2
    z_new[,,tt] <- tmp
    
  }
  return(z_new[,,-1])
  
}