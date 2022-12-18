
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


