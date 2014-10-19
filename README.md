s4vdpca
=======

```{r}
# install package using devtools
library(devtools)
install_github('mwsill/s4vdpca')

# example
library(s4vdpca)
library(broman)
p <- 1000
n <- 50
nsim <- 100
alpha <- .8
beta <- .8

Sigma <- generate_covar(alpha,beta,p)
z1 <- Sigma[[2]]
Sigma <- Sigma[[1]]

X  <- rmvn(n,rep(0,p),Sigma)
X  <- scale(X, center = TRUE, scale = FALSE)
#s4vdpca
res.s4vd1 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='bic')
res.s4vd2 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='gic2')
res.s4vd3 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='gic3')
res.s4vd4 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='gic4')
res.s4vd5 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='gic5')
res.s4vd6 <- s4vdpca(X,B=100,size=0.5,cores=4,weakness=0.5,lq=0.5,steps=500,ic_type='gic6')

plot(res.s4vd1$ic,pch='.',ylim=c(500000,800000),cex=2)
abline(v=res.s4vd1$minic)
points(res.s4vd2$ic,pch='.',col='red',cex=2)
abline(v=res.s4vd2$minic,col='red')
points(res.s4vd3$ic,pch='.',col='blue',cex=2)
abline(v=res.s4vd3$minic,col='blue')
points(res.s4vd4$ic,pch='.',col='green',cex=2)
abline(v=res.s4vd4$minic,col='green')
points(res.s4vd5$ic,pch='.',col='orange',cex=2)
abline(v=res.s4vd5$minic,col='orange')
points(res.s4vd6$ic,pch='.',col='pink')
abline(v=res.s4vd6$minic,col='pink',cex=2)
```


