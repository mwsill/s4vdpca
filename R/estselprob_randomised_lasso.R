estselprob_randomised_lasso <- function(X,u,n,lambda,B,subsets,weakness,cores){
  temp <- simplify2array(
    mclapply(1:B,mc.cores=cores,
             function(index) 
               randomised_lasso(t(X[subsets[,index],]) %*% u[subsets[,index]],lambda,weakness)!=0)
    ,higher=FALSE) 
  selprob <- rowSums(temp)/B
  avesel  <- mean(colSums(temp))  
  return(c(avesel,selprob))
}

