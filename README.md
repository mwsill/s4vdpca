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
# generate a simulated data set using the single-covariance spike model 
p <- 1000    # number of variables
n <- 50      # number of observations
alpha <- .5  # spike index 
beta <- .5   # sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha,beta,p)

# extract first eigenvector
z1 <- Sigma[[2]]

# extract variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate normal distribution using Cholesky decomposition
# see ?rmvn in package broman for details
D <- chol(Sigma)
x <- matrix(rnorm(n * p), ncol = p) %*% D + rep(rep(0,p), rep(n, p))

# apply S4VPCA with GIC5 
res <- s4vdpca(x, center=TRUE, cores=1, ic_type='gic5')

# plot the information criterion
plot(res$ic, xlab='number of features', ylab='GIC 5')
abline(v=res $minic, col='red')

# calculate angle between estimated sparse loadings vector and simulated eigenvector
angle(res$v,z1)

# calculate number of falsely selected features
sum(which(res$v!=0) %in% which(z1==0))

# apply regular pca and calculate angle
pca <- prcomp(x)
angle(pca$rotation[,1],z1)
```

