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
s4vdpca(x, center = TRUE, scale = FALSE, B = 100, size = 0.5,\cr 
        cores = 1, weakness = 0.5, lambda = NULL, nlambda = 5, \cr
        steps = 100, ic_type = "gic5")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  a numeric data matrix object, where columns represent variables and rows are observations.    
}
  \item{center}{
  Either a logical value or a numeric vector of length equal to the number of columns of \code{x}.
  The value of \code{center} determines how column centering is performed. If \code{center} is a numeric vector with length equal to the number of columns of \code{x}, then each column of \code{x} has the corresponding value from center subtracted from it. If \code{center} is \code{TRUE} then centering is done by subtracting the column means (omitting \code{NA}s) of \code{x} from their corresponding columns, and if \code{center} is \code{FALSE}, no centering is done. See generic function \code{?scale} for details.
  
}
  \item{scale}{
  Either a logical value or a numeric vector of length equal to the number of columns of \code{x}. The value of \code{scale} determines how column scaling is performed (after centering). If \code{scale} is a numeric vector with length equal to the number of columns of \code{x}, then each column of \code{x} is divided by the corresponding value from \code{scale}. If \code{scale} is \code{TRUE} then scaling is done by dividing the (centered) columns of \code{x} by their standard deviations if \code{center} is \code{TRUE}, and the root mean square otherwise. If \code{scale} is \code{FALSE}, no scaling is done. See generic function \code{?scale} for details. 
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
\eqn{\lambda_{i}^{reweighted} = \frac{\lambda}{\mathcal{U}(1-\textit{weakness},1)}}
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
Implementation of the S4VDPCA algorithm to estimate parameter estimation and feature selection consistent sparse PCs in high-dimensional, low sample size data. In a first step features are ranked by selection probabilities estimated by applying a subsampling scheme motivated by stability selection. In the second step a sparse PC is estimated by simple forward selection using one of several implemented information criteria to select an optimal model.... 

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
p <- 5000    # number of variables
n <- 50      # number of observations
alpha <- .6  # spike index 
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
set.seed(02112014)
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

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{sparse PCA}
\keyword{sparse SVD}
\keyword{stability selection}
% __ONLY ONE__ keyword per line
