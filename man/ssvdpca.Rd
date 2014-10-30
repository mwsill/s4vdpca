\name{ssvdpca}
\alias{ssvdpca}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
sparse svd pca Lee et al. 2010
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ssvdpca(X, threu = 1, threv = 1, gamu = 0, gamv = 0, u0 = svd(X)$u[, 1], v0 = svd(X)$v[, 1], merr = 10^(-4), niter = 100)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{threu}{
%%     ~~Describe \code{threu} here~~
}
  \item{threv}{
%%     ~~Describe \code{threv} here~~
}
  \item{gamu}{
%%     ~~Describe \code{gamu} here~~
}
  \item{gamv}{
%%     ~~Describe \code{gamv} here~~
}
  \item{u0}{
%%     ~~Describe \code{u0} here~~
}
  \item{v0}{
%%     ~~Describe \code{v0} here~~
}
  \item{merr}{
%%     ~~Describe \code{merr} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (X, threu = 1, threv = 1, gamu = 0, gamv = 0, u0 = svd(X)$u[, 
    1], v0 = svd(X)$v[, 1], merr = 10^(-4), niter = 100) 
{
    n = dim(X)[1]
    d = dim(X)[2]
    stop <- FALSE
    ud = 1
    vd = 1
    iter = 0
    SST = sum(X^2)
    while (vd > merr) {
        iter = iter + 1
        z = t(X) \%*\% u0
        winv = abs(z)^gamv
        sigsq = abs(SST - sum(z^2))/(n * d - d)
        tv = sort(c(0, abs(z * winv)))
        rv = sum(tv > 0)
        Bv = rep(1, d + 1) * Inf
        for (i in 1:rv) {
            lvc = tv[d + 1 - i]
            temp1 = which(winv != 0)
            temp2 = thresh(z[temp1], type = threv, delta = lvc/winv[temp1])
            vc = rep(0, d)
            vc[temp1] = temp2
            Bv[i] = sum((X - u0 \%*\% t(vc))^2)/sigsq + i * log(n * 
                d)
        }
        Iv = min(which(Bv == min(Bv)))
        temp = sort(c(0, abs(z * winv)))
        lv = temp[d + 1 - Iv]
        temp2 = thresh(z[temp1], type = threv, delta = lv/winv[temp1])
        v1 = rep(0, d)
        v1[temp1] = temp2
        v1 = v1/sqrt(sum(v1^2))
        u1 = X \%*\% v1
        u1 = u1/sqrt(sum(u1^2))
        vd = sqrt(sum((v0 - v1)^2))
        if (iter > niter) {
            stop <- TRUE
            break
        }
        u0 = u1
        v0 = v1
    }
    return(list(u = u1, v = v1, iter = iter, stop = stop, BIC = Bv))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
