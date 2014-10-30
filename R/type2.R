type2 <- function(true,est) sum(which(est==0) %in% which(true!=0))
