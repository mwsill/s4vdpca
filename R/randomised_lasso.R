randomised_lasso <- function(ols,lambda,weakness){
  lambda <- lambda/runif(length(ols),1-weakness,1)
  sign(ols) * (abs(ols) >= lambda) * (abs(ols) - lambda)
} 

