ic_rspca <- function(lambdas,p,index,adaw,ols,X,n,u0,sigsq,a,type,ic_type){
  lambda <- lambdas[p + 1 - index]
  vc <- switch(type,
               soft = softthresh(ols, delta = lambda/adaw),
               hard = hardthresh(ols, delta = lambda/adaw),
               scad = scad(ols, delta = lambda/adaw,a)
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





