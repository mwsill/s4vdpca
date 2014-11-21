s4vdpca <- function(x,center=TRUE,scale=FALSE,B=100,size=.5,cores=1,weakness=.5,a=3.7,lambda=NULL,nlambda=3,steps=100,ic_type='gic5'){
  X  <- scale(x, center, scale)
  n <- nrow(X)
  p <- ncol(X)
  #if(is.null(steps)) steps <- floor(p/100)
  svdX <- svd(X,nu=1,nv=1)
  v <- svdX$v
  u <- svdX$u
  d <- svdX$d
  ols <- as.numeric(t(X)%*%u)
  SST <- sum(X^2)
  sigsq <- abs(SST - sum(ols^2))/(n * p - p)
  subsets <- sapply(1:B,function(x){sample(1:n,n*size)})
  if(is.null(lambda)){ #search for an optimal lambda
    sels <- list()
    lqs <- seq(.1,.7,len=nlambda)
    for(i in 1:length(lqs)){
      sels[[i]] <- estselprob_randomised_lasso(X,u,n,quantile(abs(ols),lqs[i]),B,subsets,weakness,cores)[-1]
    }
    i <- which.min(abs(unlist(lapply(sels,function(x)sum(x)/p))-.5))
    lambda <- quantile(abs(ols),lqs[i])
  }
  selprobs <- estselprob_randomised_lasso(X,u,n,lambda,B,subsets,weakness,cores)[-1]
  pr <- selprobs
  selprobs <-order(selprobs,decreasing=TRUE) 
  ic <- parallel_ic(X,selprobs,p,n,sigsq,cores,steps,ic_type)
  minic <- which.min(ic)
  sv  <- subset_svd(X,selprobs[1:minic])
  out <- list(u=sv$u,v=sv$v,d=sv$d,lambda=lambda,selprobs=pr,order=selprobs,ic_type=ic_type,ic=ic,minic=minic)
  return(out)
}

