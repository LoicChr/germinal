rho_gen <- function(thetas){
  N <- (1+sqrt(1+8*length(thetas)))/2
  if (any(thetas  < 0 | thetas > pi)) stop("thetas must be between 0 and pi")
  if (!all.equal(N, as.integer(N)))  stop("number of theta parameters is wrong")

  theta_mat <- matrix(0, ncol = N, nrow = N)
  theta_mat[lower.tri(theta_mat, diag = F)] <- thetas

  B_mat <- matrix(0, nrow = nrow(theta_mat), ncol = ncol(theta_mat))
  for (i in 1:nrow(B_mat)){
    for (j in 1:i){
      if (j == 1){
        B_mat[i,j] <- cos(theta_mat[i,j])
      } else if (j >= 2 & j <= (i -1)){
        B_mat[i,j] <- cos(theta_mat[i,j])*prod(sin(theta_mat[i,1:(j-1)]))
      } else{
        B_mat[i,j] <- prod(sin(theta_mat[i,1:(j-1)]))
      }

    }
  }
  rho_mat <- B_mat %*% t(B_mat)
  return(rho_mat )
}
