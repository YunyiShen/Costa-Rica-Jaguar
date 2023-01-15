
SCR23darray <-
  function(edf, tdf){
    # edf: enconter dataframe, tdf: trap dataframe
    ### Returns 3-d array "ind x trap x occasion"
    ## stolen from https://github.com/jaroyle/scrbook/blob/master/Rpackage/scrbook/R/SCR23darray.R
    
    # Check for dups
    uq<- paste(edf[,1],edf[,2],edf[,3])
    uq<- unique(uq)
    if(any(table(uq)>1)) cat("Duplicate captures (same individual, trap, occasion) present in data set, these are not used",fill=TRUE)
    
    nind<-length(unique(edf[,1]))
    ntraps<-nrow(tdf)
    nperiods<-ncol(tdf)
    per.id<- 1:ncol(tdf)
    
    ind.id<- edf[,1]
    trap.id<- edf[,3]
    
    if( length(per.id) != length(min(per.id):max(per.id)) ){
      x<- 1:nperiods
      names(x)<-as.character(per.id)
      per.id <- x[as.character(edf[,2])]
    }
    else{
      per.id<-edf[,2]
    }
    
    y<-array(0,c(nind,ntraps, nperiods))
    
    tmp<-cbind(ind.id,trap.id,per.id)
    y[tmp]<-1
    y
  }








## rest are for plotting
plot_s <- function(individual, m_fit,...) {
  s_post <- rstan::extract(m_fit, pars = "s")$s
  s1_df <- as.data.frame(s_post[, individual, ])
  s1_sf <- st_as_sf(s1_df, coords = c("V1", "V2"))
  
  #plot(grid_sf, cex = .2, col = "grey")
  #plot(trap_sf, add = TRUE)
  #plot(trap_sf[y2d[individual, ] > 0, ], pch = 19, add = TRUE)
  plot(s1_sf, col = scales::alpha("red", .05), pch = 19,...)
}


image.scale <-
  function (z, col, x, y = NULL, size = NULL, digits = 2, labels = c("breaks", 
                                                                     "ranges"))
  {
    # sort out the location
    n <- length(col)
    usr <- par("usr")
    mx <- mean(usr[1:2]); my <- mean(usr[3:4])
    dx <- diff(usr[1:2]); dy <- diff(usr[3:4])
    if (missing(x))
      x <- mx + 1.05*dx/2	# default x to right of image
    else if (is.list(x)) {
      if (length(x$x) == 2) 
        size <- c(diff(x$x), -diff(x$y)/n)
      y <- x$y[1]
      x <- x$x[1]
    } else x <- x[1]
    if (is.null(size))
      if (is.null(y)) {
        size <- 0.618*dy/n	# default size, golden ratio
        y <- my + 0.618*dy/2	# default y to give centred scale
      } else size <- (y-my)*2/n
    if (length(size)==1)
      size <- rep(size, 2)	# default square boxes
    if (is.null(y))
      y <- my + n*size[2]/2
    # draw the image scale
    i <- seq(along = col)
    rect(x, y - i * size[2], x + size[1], y - (i - 1) * size[2], 
         col = rev(col), xpd = TRUE)
    # sort out the labels
    rng <- range(z, na.rm = TRUE)
    bks <- seq(from = rng[2], to = rng[1], length = n + 1)
    bks <- formatC(bks, format="f", digits=digits)
    labels <- match.arg(labels)
    if (labels == "breaks")
      ypts <- y - c(0, i) * size[2]
    else {
      bks <- paste(bks[-1], bks[-(n+1)], sep = " - ")
      ypts <- y - (i - 0.5) * size[2]
    }
    text(x = x + 1.2 * size[1], y = ypts, labels = bks, adj =
           ifelse(size[1]>0, 0, 1), xpd = TRUE) 
  }

SCR0density<-function (obj, nx = 30, ny = 30, Xl = NULL, Xu = NULL, Yl = NULL, 
                      Yu = NULL, scalein = 100, scaleout = 100, col="gray",ncolors = 10,whichguy=NULL) 
{
  Sout <- rstan::extract(obj, pars = "s")$s
  Sxout <- Sout[,,1]
  Syout <- Sout[,,2]
  rm(Sout)
  #z <- obj$z
  niter <- nrow(Sxout)
  if (is.null(Xl)) {
    Xl <- min(Sxout) * 0.999
    Xu <- max(Sxout) * 1.001
    Yl <- min(Syout) * 0.999
    Yu <- max(Syout) * 1.001
  }
  xg <- seq(Xl, Xu, , nx)
  yg <- seq(Yl, Yu, , ny)
  guy<-col(Sxout)
  #Sxout <- cut(Sxout[z == 1], breaks = xg)
  #Syout <- cut(Syout[z == 1], breaks = yg)
  Sxout <- cut(Sxout, breaks = xg)
  Syout <- cut(Syout, breaks = yg)
  if(is.null(whichguy)){
    Dn <- table(Sxout, Syout)/niter
    area <- (yg[2] - yg[1]) * (xg[2] - xg[1]) * scalein
    Dn <- (Dn/area) * scaleout
  }
  else{
    Dn<-table(Sxout[guy==whichguy],Syout[guy==whichguy] )/niter
  }
  
  cat("mean: ", mean(Dn), fill = TRUE)
  par(mar = c(3, 3, 3, 6))
  if (col == "gray") {
    cc <- seq(3, 17, , 10)/20
    cc <- rev(gray(cc))
  }
  else cc <- terrain.colors(ncolors)
  
  image(xg, yg, Dn, col = cc)
  image.scale(Dn, col = cc)
  box()
  return(list(grid = cbind(xg, yg), Dn = Dn))
}

