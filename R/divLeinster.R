divLeinster <- function(spxp, Z=NULL, q=2, check = TRUE){
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows and columns are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }
  spxp <- sweep(spxp, 1, rowSums(spxp), "/")
  Zp <- Z %*% t(spxp)

  if (q != 1 & q != Inf){
    mat <- t(spxp) * (Zp)^(q-1)
    mat[is.na(mat)] <- 0
    D <- colSums(mat) ^ (1/(1-q))
  }
  if (q==Inf)  {
    D <- 1/ apply(Zp, 2, max)
  }
  if (q == 1){
    D <- apply(Zp^t(spxp), 2, function(x) 1/prod(x))
  }
  return(D)
}
