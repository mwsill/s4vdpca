
screeplot.spca <- function(res,...){
  barplot(unlist(lapply(res,function(x)x$sdev))^2,ylab='Variances')
}

plot.spca <- function(res,K=1,...){
  plot(res[[K]],...)
}

plot.spc <- function(res,xlim=c(),zoom=FALSE,zf=0.01,...){
  if(zoom){
    zf <- floor(length(res$ic)*zf) 
    if(res$minic<zf) zf <- res$minic 
    plot(res$ic,ylab=paste(res$ic_type),
         ylim=c(res$ic[res$minic],max(res$ic[(res$minic-zf):(res$minic+zf)],na.rm=T))
         ,xlim=c(res$minic-zf,res$minic+zf),...)
  }else{
  plot(res$ic,ylab=paste(res$ic_type),...) 
  }
  abline(v=res$minic, col='red')
  legend("topright",paste(res$minic),text.col="red")
}