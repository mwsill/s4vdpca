#original SSVD PCA Lee et al. 2011 bzw. Shen 2006
ssvdpca <- function(X, #threu = 1, # no penalization of the left singular value
                    threv = 1,
                    #gamu = 0,  # no penalization of the left singular value
                    gamv = 0, u0 = svd(X)$u[,1], v0 = svd(X)$v[, 1], merr = 10^(-4), niter = 100)
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
    #cat("iter: ", iter, "\n")
    #cat("v: ", length(which(v0 != 0)), "\n")
    #cat("u: ", length(which(u0 != 0)), "\n")
    z = t(X) %*% u0
    winv = abs(z)^gamv
    sigsq = abs(SST - sum(z^2))/(n * d - d)
    tv = sort(c(0, abs(z * winv)))
    rv = sum(tv > 0)
    Bv = rep(1, d + 1) * Inf
    for (i in 1:rv) {
      lvc = tv[d + 1 - i]
      temp1 = which(winv != 0)
      temp2 = thresh(z[temp1], type = threv, delta = lvc/winv[temp1],a)
      vc = rep(0, d)
      vc[temp1] = temp2
      Bv[i] = sum((X - u0 %*% t(vc))^2)/sigsq + i * log(n * d)
    }
    Iv = min(which(Bv == min(Bv)))
    temp = sort(c(0, abs(z * winv)))
    lv = temp[d + 1 - Iv]
    temp2 = thresh(z[temp1], type = threv, delta = lv/winv[temp1])
    v1 = rep(0, d)
    v1[temp1] = temp2
    v1 = v1/sqrt(sum(v1^2))
    #cat("v1", length(which(v1 != 0)), "\n")
    #z = X %*% v1
    #winu = abs(z)^gamu
    #sigsq = abs(SST - sum(z^2))/(n * d - n)
    #tu = sort(c(0, abs(z * winu)))
    #ru = sum(tu > 0)
    #Bu = rep(1, n + 1) * Inf
    #for (i in 1:ru) {
    #  luc = tu[n + 1 - i]
    #  temp1 = which(winu != 0)
    #  temp2 = thresh(z[temp1], type = threu, delta = luc/winu[temp1])
    #  uc = rep(0, n)
    #  uc[temp1] = temp2
    #  Bu[i] = sum((X - uc %*% t(v1))^2)/sigsq + i * log(n * d)
    #}
    #Iu = min(which(Bu == min(Bu)))
    #temp = sort(c(0, abs(z * winu)))
    #lu = temp[n + 1 - Iu]
    #temp2 = thresh(z[temp1], delta = lu/winu[temp1])
    #u1 = rep(0, n)
    #u1[temp1] = temp2
    u1 = X %*% v1
    u1 = u1/sqrt(sum(u1^2))
    #ud = sqrt(sum((u0 - u1)^2))
    vd = sqrt(sum((v0 - v1)^2))
    if (iter > niter) {
      # print("Fail to converge! Increase the niter!")
      stop <- TRUE
      break
    }
    u0 = u1
    v0 = v1
  }
  return(list(u = u1, v = v1, d= as.numeric(t(u1)%*%X%*%v1), iter = iter, stop = stop,BIC=Bv))
}

### function to implement the soft-, hard, SCAD thresholding rule
# Input variables:
#  z: argument
#  type: thresholding rule
#    1 = (Adaptive) LASSO (default)
#		2 = hard thresholding
#		3 = SCAD
# 	delta: thresholding level
# 	a: default choice for SCAD penalty


thresh <- function(z,delta,  type=1, a=3.7){
  if(type==1){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta))
  }
  
  if(type==2){
    return(z*(abs(z)>delta))
  }
  
  if(type==3){
    return(sign(z)*(abs(z)>=delta)*(abs(z)-delta)*(abs(z)<=2*delta)+
             ((a-1)*z-sign(z)*a*delta)/(a-2)*(2*delta<abs(z))*(abs(z)<=a*delta)+z*(abs(z)>a*delta))
  }
}
