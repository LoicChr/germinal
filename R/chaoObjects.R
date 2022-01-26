
chaoObjects <- function(spxp, phy){
  if (!inherits(phy, "phylo")){
    stop("object \"phy\" is not of class \"phylo\"")}
  if (!inherits(spxp, "matrix")) {
    stop("object \"spxp\" is not of class \"matrix\"")}
  if (ncol(spxp) != length(phy$tip.label)){
    stop("object \"phy\" and object \"spxp\" does not have the same number of species")}

  Ancestors.sp <- lapply(1:length(phy$tip.label), function(x) c(Ancestors(phy, x, "all"), x))
  Branches <- lapply(Ancestors.sp, function(x) which((phy$edge[,1] %in% x) & (phy$edge[,2] %in% x)))
  Li <- unlist(lapply(Branches, function(x) sum(phy$edge.length[x]))) #Tip - root distances

  ultra <- all.equal.numeric(var(unlist(Li)), 0, tolerance = 1e-7) #Is it ultrametric
  if (!ultra) stop ("object \"phy\" must be an ultrametric tree")

  freq.dummy <- unlist(lapply(1:length(Branches),function(i) rep(i,length(Branches[[i]]))))
  desc <- lapply(unlist(Branches), function(x) Descendants(phy, phy$edge[x,2], type ="tips")[[1]])

  tmp <- lapply(desc, function(desc_i){
    x <- rep(0, length(freq.dummy))
    x[freq.dummy %in% desc_i] <- 1
    return(x)
  })
  Z <- do.call(rbind, tmp)
  pi <- sweep(spxp[,freq.dummy], 2, phy$edge.length[unlist(Branches)], FUN = "*")/Li[1]

  res <- list(Z = Z,
              pi = pi)
  return(res)
}
