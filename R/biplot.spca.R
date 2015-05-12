biplot.spca <- function(res,alpha=0.5,pcol='black',acol='gray',pcex=1,acex=1,...){
  ptx <- res[[1]]$u*(res[[1]]$d^alpha)
  pty <- res[[2]]$u*(res[[2]]$d^alpha)
  arx <- res[[1]]$v*(res[[1]]$d^(1-alpha))
  ary <- res[[2]]$v*(res[[2]]$d^(1-alpha))
  plot(ptx,pty,type="n",axes=T
       ,ylab=paste("PC2 /",sum(ary!=0),"variables",sep=" ")
       ,xlab=paste("PC1 /",sum(arx!=0),"variables",sep=" ")
       ,main='',lty=2, xlim=c(min(c(ptx,arx)),max(c(ptx,arx))),
       ylim=c(min(c(pty,ary)),max(c(pty,ary))),
       ...)
  suppressWarnings(
  arrows(0, 0, x1 = arx,y1=ary,length=0.05,col = acol)
  )
  points(ptx,pty,pch=19,cex=pcex,col=pcol)
}