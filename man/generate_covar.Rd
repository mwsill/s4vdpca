\name{generate_covar}
\alias{generate_covar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{generate a variance-covariance matrix according to single-covariance spike model described by 
Amini and Wainwright (2008)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
generate_covar(alpha, beta, p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{alpha}{
  spike index
}
  \item{beta}{
  sparsity index 
}
  \item{p}{
  number of variables
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list of
\item{Sigma}{a simulated variance-covariance matrix}
\item{z1}{the first eignvector of Sigma}
}
\references{
Sill, M., Saadati, M., and Benner, A. (2015). Applying stability selection to consistently estimate sparse principal components in high-dimensional molecular data. \cr
Oxford Bioinformatics. \cr

Amini, A.A. and Wainwright, M.J. (2008). High-dimensional analysis of semi-definite relaxations for sparse principal components.\cr
Information Theory, IEEE International Symposium, 2454-2458.  
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
n <- 100      # number of observations
alpha <- 0.8  # spike index 
beta <- 0.8   # sparsity index 

# generate a population variance covariance matrix
Sigma <- generate_covar(alpha, beta, p)

# extract first eigenvector
z1 <- Sigma[[2]]

# extract variance covariance matrix
Sigma <- Sigma[[1]]

# sample from multivariate normal distribution using Cholesky decomposition
# see ?rmvn in package broman for details
D <- chol(Sigma)
set.seed(12345)
x <- matrix(rnorm(n * p), ncol = p) %*% D + rep(rep(0,p), rep(n, p))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
