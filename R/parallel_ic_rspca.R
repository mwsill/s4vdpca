#calculate BIC stepwise in parallel break at first minimum
parallel_ic_rspca <- function(X,lambdas,p,adaw,ols,sigsq,n,u0,cores,steps,a,type,ic_type){
  ic = rep(NA, p + 1)
  points <- floor(seq(1,p,length.out=steps))
  ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
    ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
  }))
  for(i in 1:20){
    minip <- which.min(ic) 
    minp <- which.min(abs(points-which.min(ic)))
    # cat(abs(max(1,points[minp-1])-min(points[minp+1],p)), "\n")
    if(abs(max(1,points[minp-1])-min(points[minp+1],p,na.rm=T)) < steps){
      points <- max(1,points[minp-1]):min(points[minp+1],p,na.rm=T)
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
      }))
      break
    }else{
      points <- floor(seq(max(1,points[minp-1]),min(points[minp+1],p,na.rm=T),length.out=steps))
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
      }))
    }
  }
  return(ic)
}

