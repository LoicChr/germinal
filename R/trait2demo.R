
Trait2Demo <- function(tr, phis, a,b, exp = F){
  if (ncol(tr) != length(phis)+1)  stop("tr and phis don't match")
  
  xis <- numeric(length = ncol(tr))
  
  xis[1] <- cos(phis[1])
  for (i in 2:(ncol(tr)-1)) xis[i] <- prod(sin(phis[1:(i-1)]))*cos(phis[i])
  xis[ncol(tr)] <- prod(sin(phis))
  
  
  if (exp){
    theo <- exp(a*(tr %*% as.matrix(xis))+b)
  }else{
    theo <- a*(tr %*% as.matrix(xis))+b
  }
  return(as.numeric(theo))
}

