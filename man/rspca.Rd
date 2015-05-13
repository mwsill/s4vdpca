\name{rspca}
\alias{rspca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Regularized sparse prinicipal component analysis
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
rspca(x, gamv = 0, type = "soft", ic_type = "gic5" 
      , a = 3.7, merr = 10^(-4), niter = 100, cores = 1, steps = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a numeric data matrix object, where columns represent variables and rows are observations.    
}
  \item{gamv}{
weight parameter in adaptive lasso for the  
              right singular vector, nonnegative constant (default = 0, for lasso)
}
  \item{type}{
"soft" for soft thresholding, i.e. lasso and adaptive lasso \cr
"scad" for the smoothly clipped absolute deviation (SCAD) penalty \cr
"hard" for hard thresholding \cr
}
  \item{ic_type}{
Type of information criterion used for model selecion. \code{"bic"} will use the bayesian information criterion (BIC). 
Choosing \code{"gic2"} to \code{"gic6"} will apply one of the generalized information criteria (GIC) according to Yongdai et al. 2012.
}
  \item{a}{
tuning parameter for the SCAD penalty (default = 3.7)
}
  \item{merr}{
threshold to decide convergence (default = 10^(-4))
}
  \item{niter}{
maximum number of iterations (default = 100)
}
  \item{cores}{
Number of cores used for parallelization. Currently only Unix operating systems are supported.
}
  \item{steps}{
Number of points at which the information criterion is evaluted in each iteration of the parallelized optimzation algorithm to search for a global \code{ic} minimum. 
}
}
\details{
The function implements the regularized sparse PCA method proposed by Shen and Huang (2008) and Lee et al. (2010). The algorithm has been optimized for computational efficiency, i.e. to find the minimal BIC a grid search is performed instead of evaluating the BIC for each variable. In addition, this grid search algorithm hase been parallelized. 
Alternatively, different generalized information criteria (GIC) proposed to by Yongdai et al. (2012) can be used.

}
\value{
A list of
\item{u}{left singular vector}
\item{v}{sparse right singular vector, i.e. sparse loadings vector scaled to unit length}
\item{d}{singular value}
\item{iter}{number of iterations until convergence}
\item{ic_type}{type of information criterion used for model selection}
\item{ic}{vector of length \code{p}. Information criterion calculated during forward selection.}
\item{minic}{number of selected features}
}
\references{
Lee, M., Shen, H., Huang, J. Z., and Marron, J. S. (2010). Biclustering via sparse singular value decomposition.\cr
Biometrics, 66(4), 1087–1095.\cr

Shen, H. and Huang, J. (2008). Sparse principal component analysis via regularized low rank matrix approximation.\cr
Journal of Multivariate Analysis, 99(6), 1015–1034. \cr

Yongdai, K., Sunghoon, K., and Hosik, C. (2012). Consistent model selection criteria on high dimensions.\cr 
Journal of Machine Learning Research, 13, 1037–1057.\cr
}
\author{
Martin Sill m.sill\code{(at)}dkfz.de
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# generate a simulated data set using the single-covariance spike model 
p <- 1000    # number of variables
n <- 20      # number of observations
alpha <- 0.9  # spike index 
beta <- 0.7   # sparsity index 

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
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sparse PCA}
\keyword{sparse SVD}