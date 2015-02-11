aveloadings_randomised_lasso <- function(X,u,n,lambda,B,subsets,weakness,cores){
  temp <- simplify2array(
    mclapply(1:B,mc.cores=cores,
             function(index) 
               randomised_lasso(t(X[subsets[,index],]) %*% u[subsets[,index]],lambda,weakness))
    ,higher=FALSE) 
  aveloadings <- rowSums(temp)/B
  avesel  <- mean(colSums(temp!=0))  
  return(c(avesel,aveloadings))
}

