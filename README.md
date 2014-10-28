s4vdpca
=======

This R package implements methods for sparse principal component analysis as descriped in the manuscript
"Applying Stability Selection to Consistently Estimate Sparse Principal Components in High-Dimensional Molecular Data" submitted to Oxford Bioinformatics.

```{r}
# install package using devtools
library(devtools)
install_github('mwsill/s4vdpca')

# example
library(s4vdpca)
# library broman implements rmvn to simulate from a multivariate normal distribution
library(broman)

# generate a simulated data set using the single-covariance spike model 
p <- 1000
n <- 50
nsim <- 100
alpha <- .8  #spike index 
beta <- .8   #sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha,beta,p)
# extract first eigenvector
z1 <- Sigma[[2]]
# variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate data
X  <- rmvn(n,rep(0,p),Sigma)

#apply s4vdpca using different information criteria bic to gic6
res.s4vd1 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='bic')
res.s4vd2 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='gic2')
res.s4vd3 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='gic3')
res.s4vd4 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='gic4')
res.s4vd5 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='gic5')
res.s4vd6 <- s4vdpca(X,B=100,center=TRUE,size=0.5,cores=1,weakness=0.5,lq=0.5,steps=500,ic_type='gic6')

plot(res.s4vd1$ic,pch='.',ylim=c(50000,80000),xlim=c(0,400),cex=2,ylab='ic')
abline(v=res.s4vd1$minic)
points(res.s4vd1$ic,pch='.',col='black',cex=5)
points(res.s4vd2$ic,pch='.',col='red',cex=5)
abline(v=res.s4vd2$minic,col='red')
points(res.s4vd3$ic,pch='.',col='blue',cex=5)
abline(v=res.s4vd3$minic,col='blue')
points(res.s4vd4$ic,pch='.',col='green',cex=5)
abline(v=res.s4vd4$minic,col='green')
points(res.s4vd5$ic,pch='.',col='orange',cex=5)
abline(v=res.s4vd5$minic,col='orange')
points(res.s4vd6$ic,pch='.',col='pink',cex=5)
abline(v=res.s4vd6$minic,col='pink')
legend('topright',c('BIC','GIC2','GIC3','GIC4','GIC5','GIC6'),fill=c('black','red','blue','green','orange','pink'))
```

