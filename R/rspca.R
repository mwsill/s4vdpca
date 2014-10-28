#SSVD PCA Lee et al. 2011 Shen 2006
rspca <- function(X, center=TRUE, scale=FALSE, gamv = 0, type='soft', ic_type='gic5', a=3.7,  merr = 10^(-4), niter=100,cores=1,steps=100) 
{
  X  <- scale(X, center, scale)
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
    #cat("iter",iter, "\n")
    ols <- t(X) %*% u0
    adaw <- abs(ols)^gamv #adaptive lasso weights
    sigsq = abs(SST - sum(ols^2))/(n * p - p)
    lambdas = sort(c(0, abs(ols * adaw)))
    ic <- parallel_ic_rspca(X,lambdas,p,adaw,ols,sigsq,n,u0,cores,steps,a,type,ic_type)
    minbic <- which.min(ic)
    lambda <- lambda <- lambdas[p + 1 - minbic]
    v1 <- switch(type,
                 soft = softthresh(ols, delta = lambda/adaw),
                 hard = hardthresh(ols, delta = lambda/adaw),
                 scad = hardthresh(ols, delta = lambda/adaw,a)
    )     
    v1 = v1/sqrt(sum(v1^2))
    u1 = X %*% v1
    u1 = u1/sqrt(sum(u1^2))
    #ud = sqrt(sum((u0 - u1)^2))
    vd = sqrt(sum((v0 - v1)^2))
    if (iter > niter) {
      # print("Fail to converge! Increase the niter!")
      stop <- TRUE
      break
    }
    u0 = u1
    v0 = v1
  }
  return(list(u = u1, v = v1, d= as.numeric(t(u1)%*%X%*%v1), iter = iter, ic_type=ic_type, ic=ic , stop = stop))
}
