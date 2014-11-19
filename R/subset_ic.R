subset_ic <- function(X,subset,n,p,sigsq,ic_type){
  sv  <- subset_svd(X,subset)
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

