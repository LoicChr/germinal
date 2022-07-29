findLVeq <- function(K, A){
  alphas.inv <- pseudoinverse(A)
  Ns <- alphas.inv %*% as.matrix(K)
  # at least one species has been competitively excluded
  excluded <- c()
  if(any(Ns < 0)){
    excluded <- c(excluded, which.min(Ns))
    Ks_red <- K[-excluded]
    alphas <- A[-excluded, -excluded, drop = F]
    Ns <- K*0
    Ns[-excluded] <- pseudoinverse(alphas) %*% as.matrix(Ks_red)

    while(any(Ns < 0)){
      excluded <- c(excluded, which.min(Ns))
      Ks_red <- K[-excluded]
      alphas <- as.matrix(A[-excluded, -excluded, drop = F])
      Ns <- K*0
      Ns[-excluded] <- pseudoinverse(alphas) %*% as.matrix(Ks_red)
    }
  }
  return(Ns)
}
