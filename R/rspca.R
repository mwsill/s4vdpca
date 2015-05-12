#SSVD PCA Lee et al. 2011 Shen 2006
rspca <- function(X, gamv = 0, type='soft', ic_type='gic5', a=3.7,  merr = 10^(-4), niter=100,cores=1,steps=100) 
{
  n <- nrow(X)
  p <- ncol(X)
  stop <- FALSE
  svdX <- svd(X,nu=1,nv=1)
  v0 <- svdX$v
  u0 <- svdX$u
  vd <- 1
  iter <- 0
  SST <- sum(X^2)
  while (vd > merr) {
    iter <- iter + 1
    ols <- t(X) %*% u0
    adaw <- abs(ols)^gamv #adaptive lasso weights
    sigsq = abs(SST - sum(ols^2))/(n * p - p)
    lambdas = sort(c(0, abs(ols * adaw)))
    ic <- .parallel_ic_rspca(X,lambdas,p,adaw,ols,sigsq,n,u0,cores,steps,a,type,ic_type)
    minbic <- which.min(ic)
    lambda <- lambda <- lambdas[p + 1 - minbic]
    v1 <- switch(type,
                 soft = .softthresh(ols, delta = lambda/adaw),
                 hard = .hardthresh(ols, delta = lambda/adaw),
                 scad = .scad(ols, delta = lambda/adaw,a)
    )     
    v1 = v1/sqrt(sum(v1^2))
    u1 = X %*% v1
    u1 = u1/sqrt(sum(u1^2))
    #ud = sqrt(sum((u0 - u1)^2))
    vd = sqrt(sum((v0 - v1)^2))
    if (iter > niter) {
      print("Fail to converge! Increase the niter!")
      stop <- TRUE
      break
    }
    u0 = u1
    v0 = v1
  }
  minic <- which.min(ic)
  d <- as.numeric(t(u1)%*%X%*%v1)
  out <- list(u = u1, v = v1, d=d ,sdev=d/sqrt(n-1), iter = iter, ic_type=ic_type, ic=ic, minic=minic)
  class(out) <- 'spc'
  structure(out,call=match.call())
}

.softthresh <- function(ols,delta) return(sign(ols) * (abs(ols) >= delta) * (abs(ols) - delta))

.scad <- function(ols, delta, a){ 
  return(sign(ols) * (abs(ols) >= delta) * (abs(ols) - delta) * 
           (abs(ols) <= 2 * delta) + ((a - 1) * ols - sign(ols) * 
                                        a * delta)/(a - 2) * (2 * delta < abs(ols)) * (abs(ols) <= 
                                                                                         a * delta) + ols * (abs(ols) > a * delta)) 
}

.hardthresh <- function(ols,delta) return(ols * (abs(ols) > delta))

.ic_rspca <- function(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type){
  lambda <- lambdas[p + 1 - index]
  vc <- switch(type,
               soft = .softthresh(ols, delta = lambda/adaw),
               hard = .hardthresh(ols, delta = lambda/adaw),
               scad = .scad(ols, delta = lambda/adaw,a)
  )           
  ic <- switch(ic_type,
               bic = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 log(p*n),
               gic2 = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 p^(1/3),
               gic3 = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 2*log(p),
               gic4 = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 2*log(p+log(log(p))),
               gic5 = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 log(log(n*p))*log(p), 
               gic6 = sum((X - u0 %*% t(vc))^2)/sigsq + index *
                 log(n*p)*log(p) 
  )
}

.parallel_ic_rspca <- function(X,lambdas,p,adaw,ols,sigsq,n,u0,cores,steps,a,type,ic_type){
  ic = rep(NA, p + 1)
  points <- floor(seq(1,p,length.out=steps))
  ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
    .ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
  }))
  for(i in 1:20){
    minip <- which.min(ic) 
    minp <- which.min(abs(points-which.min(ic)))
    if(abs(max(1,points[minp-1])-min(points[minp+1],p,na.rm=T)) < steps){
      points <- max(1,points[minp-1]):min(points[minp+1],p,na.rm=T)
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        .ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
      }))
      break
    }else{
      points <- floor(seq(max(1,points[minp-1]),min(points[minp+1],p,na.rm=T),length.out=steps))
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        .ic_rspca(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type)
      }))
    }
  }
  return(ic)
}
