scad <- function(ols, delta, a){ 
  return(sign(ols) * (abs(ols) >= delta) * (abs(ols) - delta) * 
           (abs(ols) <= 2 * delta) + ((a - 1) * ols - sign(ols) * 
                                        a * delta)/(a - 2) * (2 * delta < abs(ols)) * (abs(ols) <= 
                                                                                         a * delta) + ols * (abs(ols) > a * delta)) 
}