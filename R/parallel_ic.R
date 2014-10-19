# BIC shortcut
parallel_ic <- function(X,selprobs,p,n,sigsq,cores,steps,ic_type){
  ic <- rep(NA,p)
  points <- floor(seq(1,p,length.out=steps))
  ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
    subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
  }))
  for(i in 1:20){
    minip <- which.min(ic) 
    minp <- which.min(abs(points-which.min(ic)))
    if(abs(max(1,points[minp-1])-min(points[minp+1],p,na.rm=T)) < steps){
      points <- max(1,points[minp-1]):min(points[minp+1],p,na.rm=T)
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
      }))
      break
    }else{
      points <- floor(seq(max(1,points[minp-1]),min(points[minp+1],p,na.rm=T),length.out=steps))
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
      }))
    }
  }
  return(ic)
}

