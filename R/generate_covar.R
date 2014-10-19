#generate Sigma covaiance spike model
generate_covar <- function(alpha,beta,p){
  # first eigen value
  lambda1 <- p^alpha  #first eigenvalue std
  # first eigenvector
  p1 <- floor(p^beta)
  # first eigenvector
  z1 <- c(rep(1/sqrt(p1),p1),rep(0,p-p1))
  Sigma <- (lambda1-1) * (z1%*%t(z1)) +diag(p)
  return(list(Sigma,z1))
}  