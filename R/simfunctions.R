# functions to perform simulation study
#
# Martin Sill
# m.sill@dkfz.de
# 15.08.2013

genSigmak2 <- function(alpha,beta,gamma,p){
  lambda1 <- p^alpha  #first eigenvalue std
  lambda2 <- p^(alpha-gamma) 
  p1 <- floor(p^beta)
  p2 <- floor((p-p1)^(beta-gamma))            
  z1 <- c(rep(1/sqrt(p1),p1),rep(0,p-p1))
  # second
  z2 <- c(rep(0,p1),rep(1/sqrt(p2),p2),rep(0,p-p1-p2))
  # third
  #z3 <- c(rep(0,p1+p2),rep(1/sqrt(p3),p3),rep(0,p-p1-p2-p3))
  Sigma <- (lambda1-1) * (z1%*%t(z1)) + (lambda2-1) * (z2%*%t(z2)) + diag(p)
  return(list(Sigma,z1,z2))
}  

#generate Sigma covaiance spike model
genSigma <- function(alpha,beta,p){
  # first eigen value
  lambda1 <- p^alpha  #first eigenvalue std
  # first eigenvector
  p1 <- floor(p^beta)
  # first eigenvector
  z1 <- c(rep(1/sqrt(p1),p1),rep(0,p-p1))
  Sigma <- (lambda1-1) * (z1%*%t(z1)) +diag(p)
  return(list(Sigma,z1))
}  

#calculate angle 
angle <- function(a,b,digits=10){
  dotproduct <- round(abs(sum(a*b)),digits=10)
  acos(dotproduct) * (180/pi)
}

# calculate type 1 errors
type1 <- function(true,est) sum(which(est!=0) %in% which(true==0))

# calculate type 2 errors
type2 <- function(true,est) sum(which(est==0) %in% which(true!=0))

imagett <- function(x,...){
  x[which(is.nan(x),arr.ind=T)] <- 0
  imagett(x,...)
}
