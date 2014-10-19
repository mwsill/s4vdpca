softthresh <- function(ols,delta) return(sign(ols) * (abs(ols) >= delta) * (abs(ols) - delta))


