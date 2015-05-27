s4vdpca
=======

This R package implements methods for sparse principal component analysis as described in the Bioinformatics article:

[Applying Stability Selection to Consistently Estimate Sparse Principal Components in High-Dimensional Molecular Data](http://bioinformatics.oxfordjournals.org/content/early/2015/04/28/bioinformatics.btv197.long) 

```{r}
# install package using devtools
# install.packages('devtools')
library(devtools)                  
install_github('mwsill/s4vdpca')
library(s4vdpca)

# generate a simulated data set using the single-covariance spike model 
p <- 5000    # number of variables
n <- 100      # number of observations
alpha <- 0.8  # spike index 
beta <- 0.8   # sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha,beta,p)

# extract first eigenvector
z1 <- Sigma[[2]]

# extract variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate normal distribution using Cholesky decomposition
# see ?rmvn in package broman for details
D <- chol(Sigma)
set.seed(24022015)
x <- matrix(rnorm(n * p), ncol = p) %*% D + rep(rep(0,p), rep(n, p))

#show documentation
?s4vdpca
?rspca

# mean centering
x <- scale(x,center=TRUE,scale=FALSE)

# apply S4VDPCA and RSPCA with different penalization functions, all with GIC5 
# parallelization is not yet available on Windows machines
res1 <- s4vdpca(x, cores=1, ic_type='gic5') #s4vdpca
res2 <- rspca(x, cores=1, ic_type='gic5') #lasso
res3 <- rspca(x, cores=1, ic_type='gic5', type='scad') #scad 
res4 <- rspca(x, cores=1, ic_type='gic5', gamv=1) # adaptive lasso

# plot the information criterion
par(mfrow=c(2,2))
plot(res1,main='S4VDPCA')
plot(res2,main='RSPCA lasso')
plot(res3,main='RSPCA scad')
plot(res4,main='RSPCA adaptive lasso')
```
![](https://github.com/mwsill/s4vdpca/blob/figures/img1.png)

```{r}
# calculate angle between estimated sparse loadings vector and simulated eigenvector
angle(res1$v,z1)
angle(res2$v,z1)
angle(res3$v,z1)
angle(res4$v,z1)

# calculate number of falsely selected features
type1(z1,res1$v)
type1(z1,res2$v)
type1(z1,res3$v)
type1(z1,res4$v)

# calculate number of type2 errors
type2(z1,res1$v)
type2(z1,res2$v)
type2(z1,res3$v)
type2(z1,res4$v)

# calculate empirical spike index
emp_spike(res1$sdev,p)
emp_spike(res2$sdev,p)
emp_spike(res3$sdev,p)
emp_spike(res4$sdev,p)
alpha

# calculate empirical sparsity index
emp_sparsity(res1$minic,p)
emp_sparsity(res2$minic,p)
emp_sparsity(res3$minic,p)
emp_sparsity(res4$minic,p)
beta
```

Real data application 
```{r}
# loading the medulloblastoma example data set
# requires the Bioconductor package Biobase
# http://bioconductor.org/
# source("http://bioconductor.org/biocLite.R")
# biocLite()

data(medullo)

# data is stored as ExpressionSet
medullo

# extract gene expression matrix
X <- t(exprs(medullo))
dim(X)

# perform sparse PCA by S4VPCA 
res <- spca(X,K=3,center=TRUE,method='s4vdpca')

# plot results
par(mfrow=c(2,2))
screeplot(res)
plot(res,1,main='PC1')
plot(res,2,main='PC2')
biplot(res,pcol=medullo$subtype,alpha=0.4)
```
![](https://github.com/mwsill/s4vdpca/blob/figures/img2.png)

```{r}
# interactive 3d scatterplot using the threejs package
# and the wesanderson color palette generator
library(threejs)
library(wesanderson)
cols <- wes_palette("Royal1")
names(cols) <- unique(medullo$subtype)
cols <- cols[medullo$subtype]
names(cols) <- NULL
scatterplot3js(cbind(PC1=res[[1]]$u,PC2=res[[2]]$u,PC3=res[[3]]$u),
                col=cols,
                size=3,
                renderer='canvas')
```
[![](https://github.com/mwsill/s4vdpca/blob/figures/scatter3djs.png)](http://htmlpreview.github.io/?https://github.com/mwsill/s4vdpca/blob/figures/scatter3djs.html)