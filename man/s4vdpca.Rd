\name{s4vdpca}
\alias{s4vdpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse principal component analysis using stability selection.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
s4vdpca(X , B = 100, size = 0.5,
        cores = 1, weakness = 0.5, lambda = NULL, nlambda = 5,
        steps = 100, ic_type = "gic5")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
  a numeric data matrix object, where columns represent variables and rows are observations.    
}
  \item{B}{
Number of subsamples used for stability selection.
}
  \item{size}{
Proportion of samples drawn without replacement for each subsample.
}
  \item{cores}{
Number of cores used for parallelization. Currently only Unix operating systems are supported.
}
  \item{weakness}{
Parameter for the randomised lasso used within the stability selection, i.e. will reweight the penalization parameter \code{lambda} applied to each coefficient using random weights sampled from an uniform distribution. \cr
\eqn{\lambda_{i}^{reweighted} = \frac{\lambda}{\mathcal{U}(\textit{weakness},1)}}
}
  \item{lambda}{
Penalization parameter used for the pointwise stability selection.
}
  \item{nlambda}{
Number of lambdas evaluated to find optimal selection probabilities.
}
  \item{steps}{
Number of points at which the information criterion is evalute in each iteration of the parallelized optimzation algorithm to search for a global \code{ic} minimum. 
}
  \item{ic_type}{
Type of information criterion used for model selecion. \code{"bic"} will use the bayesian information criterion (BIC). 
Choosing \code{"gic2"} to \code{"gic6"} will apply one of the generalized information criteria (GIC) according to Yongdai et al. 2012.
}
}
\details{
Implementation of the S4VDPCA algorithm to estimate parameter estimation and feature selection consistent sparse PCs in high-dimensional, low sample size data. In a first step features are ranked by selection probabilities estimated by applying a subsampling scheme motivated by stability selection. In the second step a sparse PC is estimated by simple forward selection using one of several implemented information criteria to select an optimal model...

}
\value{
A list of
\item{u}{left singular vector}
\item{v}{sparse right singular vector, i.e. sparse loadings vector scaled to unit length}
\item{d}{singular value}
\item{lambda}{penalization parameter used for the stability selection}
\item{selprobs}{estimated selection probabilities}
\item{order}{ordering of the features according to the selection probabilities}
\item{ic_type}{type of information criterion used for model selection}
\item{ic}{vector of length \code{p}. Information criterion calculated during forward selection.}
\item{minic}{number of selected features}
}
\references{
Yongdai, K., Sunghoon, K., and Hosik, C. (2012). Consistent model selection criteria on high dimensions.\cr 
Journal of Machine Learning Research, 13, 1037–1057.\cr

Meinshausen N. and B\"uhlmann P. (2010). Stability Selection. \cr
Journal of the Royal Statistical Society: Series B (Statistical Methodology) 72, 4, 417-473.\cr
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
\keyword{stability selection}
% __ONLY ONE__ keyword per line
