
#calculate angle 
angle <- function(a,b,digits=10){
  dotproduct <- round(abs(sum(a*b)),digits=10)
  acos(dotproduct) * (180/pi)
}

# calculate type 1 errors
type1 <- function(true,est) sum(which(est!=0) %in% which(true==0))

# calculate type 2 errors
type2 <- function(true,est) sum(which(est==0) %in% which(true!=0))

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

emp_spike <- function(sdev,p){
  log(sdev^2)/log(p)
}

emp_sparsity <- function(nsel,p){
  log(nsel)/log(p)
}
