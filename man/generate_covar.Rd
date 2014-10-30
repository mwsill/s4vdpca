\name{generate_covar}
\alias{generate_covar}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{generate variance covariance matrix according to the spike model of ...
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
%%     ~~Describe \code{alpha} here~~
}
  \item{beta}{
%%     ~~Describe \code{beta} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
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
function (alpha, beta, p) 
{
    lambda1 <- p^alpha
    p1 <- floor(p^beta)
    z1 <- c(rep(1/sqrt(p1), p1), rep(0, p - p1))
    Sigma <- (lambda1 - 1) * (z1 \%*\% t(z1)) + diag(p)
    return(list(Sigma, z1))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
