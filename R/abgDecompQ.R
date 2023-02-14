
abgDecompQ <- function(spxp, Z=diag(ncol(spxp)), q=2, site.weight = NULL, check=TRUE) {
  #Calcul the diversity of each site of sites by species matrix.
  #spxp columns and Z rows/cols are assumed to be in the same order.

  if (check){
    if (!inherits(spxp, "matrix")) {
      stop("object \"spxp\" is not of class \"matrix\"")}
    if (!inherits(Z, "matrix")) {
      stop("object \"Z\" is not of class \"matrix\"")}
    if (!all(c(ncol(Z), nrow(Z)) == ncol(spxp))){
      stop("object \"Z\" and object \"spxp\" do not have matching dimensions")}
    if (nrow(spxp) != length(site.weight) & !is.null(site.weight)){
      stop("object \"site.weight\" and object \"spxp\" do not have matching dimensions")
    }
  }
  if (is.null(site.weight)){
    site.weight <- rep(1/nrow(spxp), nrow(spxp))
  }else{
    site.weight <- site.weight/sum(site.weight)
  }

  spxp <- sweep(spxp, 1, rowSums(spxp), "/")

  gamma.ab <- colSums(sweep(spxp, 1, site.weight, "*"))

  Gamma <- divLeinster(t(as.matrix(gamma.ab)), Z=Z , q=q, check = FALSE)
  Alphas <- divLeinster(spxp, Z=Z , q=q, check = FALSE)

  if (q != 1 & q != Inf) {
    mAlpha <- (sum(site.weight * (Alphas ^ (1 - q))))^(1 / (1 - q))
  }
  if (q==1){
    mAlpha <- exp(sum(site.weight * log(Alphas)))
  }
  if (q==Inf){
    mAlpha <- min(Alphas)
  }
  Beta <- Gamma / mAlpha

  names(Alphas) <- row.names(spxp)
  res <- list(Gamma=Gamma, Beta=Beta, mAlpha=mAlpha, Alphas=Alphas)

  return(res)
}
