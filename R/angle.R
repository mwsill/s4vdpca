angle <- function(a,b,digits=10){
  dotproduct <- round(abs(sum(a*b)),digits=10)
  acos(dotproduct) * (180/pi)
}
