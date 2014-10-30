s4vdpca
=======

This R package implements methods for sparse principal component analysis as descriped in the manuscript
"Applying Stability Selection to Consistently Estimate Sparse Principal Components in High-Dimensional Molecular Data" submitted to Oxford Bioinformatics.

```{r}
# generate a simulated data set using the single-covariance spike model 
p <- 1000    # number of variables
n <- 50      # number of observations
alpha <- .9  # spike index 
beta <- .9   # sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha,beta,p)

# extract first eigenvector
z1 <- Sigma[[2]]

# extract variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate normal distribution using Cholesky decomposition
# see ?rmvn in package broman for details
D <- chol(Sigma)
set.seed(30102014)
x <- matrix(rnorm(n * p), ncol = p) %*% D + rep(rep(0,p), rep(n, p))

# apply S4VDPCA and RSPCA with different penalization functions, all with GIC5 
res1 <- s4vdpca(x, center=TRUE, cores=4, ic_type='gic5')
res2 <- rspca(x, center=TRUE, cores=4, ic_type='gic5') #lasso
res3 <- rspca(x, center=TRUE, cores=4, ic_type='gic5', type='scad') #scad 
res4 <- rspca(x, center=TRUE, cores=4, ic_type='gic5', gamv=1) # adaptive lasso

# plot the information criterion
par(mfrow=c(2,2))
plot(res1$ic, xlab='number of selected features', ylab='GIC 5',main='S4VDPCA')
abline(v=res1$minic, col='red')
text(y=max(res1$ic,na.rm=T)-1000,x=res1$minic+100,res1$minic,col='red')
plot(res2$ic, xlab='number of selected features', ylab='GIC 5',main='RSPCA lasso')
abline(v=res2$minic, col='red')
text(y=max(res2$ic,na.rm=T)-1000,x=res2$minic+100,res2$minic,col='red')
plot(res3$ic, xlab='number of selected features', ylab='GIC 5',main='RSPCA scad')
abline(v=res3$minic, col='red')
text(y=max(res3$ic,na.rm=T)-1000,x=res3$minic+100,res3$minic,col='red')
plot(res4$ic, xlab='number of selected features', ylab='GIC 5',main='RSPCA adaptive lasso')
abline(v=res4$minic, col='red')
text(y=max(res4$ic,na.rm=T)-1000,x=res4$minic+100,res4$minic,col='red')

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

# calculate type 2 errors
type2(z1,res1$v)
type2(z1,res2$v)
type2(z1,res3$v)
type2(z1,res4$v)

# apply regular PCA and calculate angle between loadings vector
# and simulated eigenvector
pca <- prcomp(x)
angle(pca$rotation[,1],z1)

# ssvdpca is the original rspca function by Lee et al. 2010
X  <- scale(x, center=TRUE)
res5 <- ssvdpca(X) #lasso
# optimized code bic evaluted for each variable 
res6 <- rspca(X, center=FALSE, cores=1,steps=1000, ic_type='bic') #lasso
# estimated loadings are the same 
all(res5$v==res6$v)
```

