###############################################################################
#                                                                             #
#   Functions to generate the interaction matrix                               #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin       #
#   *contact                                                                  #
###############################################################################

# rho_gen calculates a correlation matrix based on angle parameters (thetas) that determines the Cholesky factorisation of the final correlation matrix.
# It was coded directly from Forrester et Zhang 2020 (Journal of Multivariate Analysis)

# Interaction_matrix calculates the pairwise competition matrix as a function of the intercept, mu, sigma and rhos parameters. Intra specifies the intraspecific interaction coefficients.


interactionMatrix <- function(tr, intercept, mu, sigma, rho = NULL, intra = 1, std = F){
  if (length(unique(c(ncol(tr), length(mu), length(sigma)))) > 1)  stop("mu, sigma and ncol(tr) don't match")
  if (length(intercept) > 1) stop("intercept must be a scalar")
  if (!(length(intercept) == nrow(tr) | length(intra) == 1 | is.null(intra))) stop("intra must be either NULL, a scalar or a vector of length equal to nrow(tr)")

  if (is.null(rho)){
    rho_mat <- diag(ncol(tr))
  }else{
    N <- (1+sqrt(1+8*length(rho)))/2
    if (N != ncol(tr)) stop("number of theta parameters don't match the number of traits")

    rho_mat <- rho_gen(acos(rho))
  }
  if (ncol(tr) == 1){
    inv_sigma_mat <- as.matrix(1/sigma)
  }else{
    inv_sigma_mat <- solve(diag(sigma) %*% rho_mat %*% diag(sigma))
  }

  aij = matrix(NA, nrow(tr), nrow(tr))
  for (i in 1:nrow(tr)){
    for (j in 1:nrow(tr)){
      dt <- as.matrix(tr[i,]- tr[j,])
      aij[i,j] = exp(-0.5* t(dt -mu) %*% inv_sigma_mat %*% (dt -mu))
    }
  }
  aij = intercept*aij

  if (std) aij = aij*exp(0.5 *t(mu) %*% inv_sigma_mat %*% (mu))

  if (!is.null(intra)) diag(aij) <- intra

  return(aij)
}
