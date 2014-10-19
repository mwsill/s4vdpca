s4vdpca <- function(X,B=100,size=.5,cores=1,weakness=.5,a=3.7,lambda=NULL,lq=0.5,steps=100,ic_type='gic5'){
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
  if(is.null(lambda)) lambda <- quantile(abs(ols),lq)
  selprobs <- estselprob_randomised_lasso(X,u,n,lambda,B,subsets,weakness,cores)[-1]
  pr <- selprobs
  selprobs <-order(selprobs,decreasing=TRUE) 
  ic <- parallel_ic(X,selprobs,p,n,sigsq,cores,steps,ic_type)
  minic <- which.min(ic)
  sv  <- subset_svd(X,selprobs[1:minic])
  out <- list(u=sv$u,v=sv$v,d=sv$d,lambda=lambda,selprobs=selprobs,ic_type=ic_type,ic=ic,minic=minic,pr=pr)
  return(out)
}

