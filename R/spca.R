spca <- function(x,K=3,center=TRUE,scale=FALSE,method='s4vdpca',...){
  if(center) x <- scale(x,center,scale) 
  out <- list()
  for(k in 1:K){
    message('estimating PC',k, appendLF = FALSE)
    switch(method,
           s4vdpca={out[[k]] <- s4vdpca(x,...)},
           rspca={out[[k]] <- rspca(x,...)}
    )
    x <- x - (out[[k]]$d * out[[k]]$u %*% t(out[[k]]$v)) 
    message(' done.') 
 }
 class(out) <- 'spca'
 attr(out, 'method') <- paste(method)
 structure(out, call=match.call())
}

summary.spca <- function(object){
  vars <- unlist(lapply(object,function(x)x$sdev))^2
  vars <- vars/sum(vars)
  importance <- cbind(`Standard deviation` = unlist(lapply(object,function(x)x$sdev))^2,
                      `Proportion of Variance` = round(vars, 2), `Cumulative Proportion` = round(cumsum(vars), 2),
                      `Non-zero Loadings` = unlist(lapply(object,function(x)x$minic)),
                      `Empirical sparsity index` = emp_sparsity(unlist(lapply(object,function(x)x$minic))
                                                                ,unlist(lapply(object,function(x)length(x$v)))),
                      `Empirical spike index` = emp_sparsity(unlist(lapply(object,function(x)x$sdev^2))
                                                             ,unlist(lapply(object,function(x)length(x$v)))),
                      `IC` = unlist(lapply(object,function(x)round(x$ic[x$minic],2))))   
  rownames(importance) <- paste('PC',1:length(object),sep='')
  colnames(importance)[7] <- paste(object[[1]]$ic_type)
  #object$importance <- importance
  #class(object) <- "summary.spca"
  cat("\nCall:\n", paste(deparse(attr(object,"call")), sep = "\n", collapse = "\n"), 
      "\n", sep = "")
  cat("\nMethod:\n", attr(object,"method"),
      "\n", sep = "")
  importance
}

print.spca <- function(object){
  print(summary.spca(object))
}