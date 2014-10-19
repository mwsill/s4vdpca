subset_svd <- function(X,subset){  
  sv  <- svd(X[,subset],nu=1,nv=1) #irlba is mostly too slow here
  v <- rep(0,ncol(X))
  v[subset] <- sv$v
  out <- list(u=sv$u,v=v,d=sv$d[1])
}