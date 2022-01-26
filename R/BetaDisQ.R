
BetaDisQ <- function(spxp, Z=NULL, q=2, check = TRUE){
  #Calcul the site pairwise diversity of a sites by species matrix.
  #spxp columns and Z rows/cols are assumed to be in the same order.
  if (is.null(Z)) Z <- diag(ncol(spxp))
  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" does not have matching dimensions")}
  }

  N <- nrow(spxp)
  dis <- matrix(NA, N, N)
  for (i in 2:N) {
    for (j in 1:(i-1)) {
      spxp.dummy <- spxp[c(i,j), ]
      res <- abgDecompQ(as.matrix(spxp.dummy), Z = Z, q = q, check = FALSE)
      dis[i, j] <- dis[j, i] <- res$Beta
    }
  }

  diag(dis) <- 1
  dis <- dis - 1
  row.names(dis) <- colnames(dis) <- row.names(spxp)
  return(dis)
}
