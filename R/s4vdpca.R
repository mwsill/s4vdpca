s4vdpca <- function(X,B=500,size=.5,cores=1,weakness=.5,a=3.7,rankbyloadings=F,lambda=NULL,nlambda=5,steps=100,ic_type='gic5'){
  n <- nrow(X)
  p <- ncol(X)
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
    lqs <- seq(.1,.9,len=nlambda)
    if(!rankbyloadings){
      for(i in 1:length(lqs)){
        sels[[i]] <- .estselprob_randomised_lasso(X,u,n,quantile(abs(ols),lqs[i]),B,subsets,weakness,cores)[-1]
      }
      i <- which.max(unlist(lapply(sels,function(x)length(unique(rank(x)))))) # minimise ties maybe TODO check hartigans dip statistic
    }else{
      for(i in 1:length(lqs)){
        sels[[i]] <- .aveloadings_randomised_lasso(X,u,n,quantile(abs(ols),lqs[i]),B,subsets,weakness,cores)[-1]
      }
      i <- which.max(unlist(lapply(sels,function(x)length(unique(rank(abs(x)))))))
    }  
    lambda <- quantile(abs(ols),lqs[i])
  }
  if(!rankbyloadings){
    selprobs <- .estselprob_randomised_lasso(X,u,n,lambda,B,subsets,weakness,cores)[-1]
  }else{
    selprobs <- .aveloadings_randomised_lasso(X,u,n,lambda,B,subsets,weakness,cores)[-1]
  }
  pr <- selprobs
  #selprobs <-order(selprobs,decreasing=TRUE) 
  selprobs <- order(rank(selprobs,ties.method='random'),decreasing=T) # rank ties at random
  if(rankbyloadings) selprobs <- order(rank(abs(pr),ties.method='random'),decreasing=T) 
  ic <- .parallel_ic(X,selprobs,p,n,sigsq,cores,steps,ic_type)
  minic <- which.min(ic)
  sv  <- .subset_svd(X,selprobs[1:minic])
  out <- list(u=sv$u,v=sv$v,d=sv$d,sdev=sv$d/sqrt(n-1),lambda=lambda,selprobs=pr,order=selprobs,ic_type=ic_type,ic=ic,minic=minic)
  class(out) <- 'spc'
  structure(out,call=match.call())
}

.randomised_lasso <- function(ols,lambda,weakness){
  lambda <- lambda/runif(length(ols),weakness-1,1)
  sign(ols) * (abs(ols) >= lambda) * (abs(ols) - lambda)
} 

.subset_ic <- function(X,subset,n,p,sigsq,ic_type){
  sv  <- .subset_svd(X,subset)
  ic <- switch(ic_type,
               bic = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 log(p+n),
               gic2 = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 p^(1/3),
               gic3 = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 2*log(p),
               gic4 = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 2*log(p+log(log(p))),
               gic5 = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 log(log(n*p))*log(p), 
               gic6 = sum((X - sv$d*sv$u%*%t(sv$v))^2)/sigsq + length(subset) *
                 log(n*p)*log(p) 
  )
  return(ic)
}

.parallel_ic <- function(X,selprobs,p,n,sigsq,cores,steps,ic_type){
  ic <- rep(NA,p)
  points <- floor(seq(1,p,length.out=steps))
  ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
    .subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
  }))
  for(i in 1:20){
    minip <- which.min(ic) 
    minp <- which.min(abs(points-which.min(ic)))
    if(abs(max(1,points[minp-1])-min(points[minp+1],p,na.rm=T)) < steps){
      points <- max(1,points[minp-1]):min(points[minp+1],p,na.rm=T)
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        .subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
      }))
      break
    }else{
      points <- floor(seq(max(1,points[minp-1]),min(points[minp+1],p,na.rm=T),length.out=steps))
      ic[points] <- unlist(mclapply(points,mc.cores=cores, function(index){
        .subset_ic(X,selprobs[1:index],n,p,sigsq,ic_type)
      }))
    }
  }
  return(ic)
}

.estselprob_randomised_lasso <- function(X,u,n,lambda,B,subsets,weakness,cores){
  temp <- simplify2array(
    mclapply(1:B,mc.cores=cores,
             function(index) 
               .randomised_lasso(t(X[subsets[,index],]) %*% u[subsets[,index]],lambda,weakness)!=0)
    ,higher=FALSE) 
  selprob <- rowSums(temp)/B
  avesel  <- mean(colSums(temp))  
  return(c(avesel,selprob))
}

.aveloadings_randomised_lasso <- function(X,u,n,lambda,B,subsets,weakness,cores){
  temp <- simplify2array(
    mclapply(1:B,mc.cores=cores,
             function(index) 
               .randomised_lasso(t(X[subsets[,index],]) %*% u[subsets[,index]],lambda,weakness))
    ,higher=FALSE) 
  aveloadings <- rowSums(temp)/B
  avesel  <- mean(colSums(temp!=0))  
  return(c(avesel,aveloadings))
}

.subset_svd <- function(X,subset){  
  sv  <- svd(X[,subset],nu=1,nv=1) #irlba is mostly too slow here
  v <- rep(0,ncol(X))
  v[subset] <- sv$v
  out <- list(u=sv$u,v=v,d=sv$d[1])
}
