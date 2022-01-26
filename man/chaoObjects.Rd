\name{chaoObjects}
\alias{chaoObjects}
\title{ Preparation of the objects necessary to compute Chao's diversity index}
\usage{
chaoObjects(spxp, phy)
}
\description{
chaoObjects is a data preparation function. chaoObjects returns adequate arguments for abgDecompQ, BetaDisQ and divLeinster to perform a diversity analysis using Chao's diversity index.

Warning: this formula only works with an ultrametric tree!
}
\arguments{
\item{spxp}{Site by species matrix}
\item{phy}{Ultrametric phylogenetic tree}
}
\references{
Chalmandrier et al. (2015) Effects of species' similarity and dominance on the functional and phylogenetic structure of a plant meta‐community. Ecology 96, 143-153

Chao et al. 2010 Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B: Biological Sciences 365: 3599-3609

Leinster & Cobbold (2012) Measuring diversity: the importance of species similarity. Ecology 93, 477-489 \\
}
\examples{
library(ape)
data(phylocom, package = "picante")
spxp <- phylocom$sample
phy <- phylocom$phylo
is.missing<- !(phy$tip.label \%in\% colnames(spxp))
phy <- drop.tip(phy, phy$tip.label[is.missing])
spxp <- spxp[, phy$tip.label]
chao.objs <- chaoObjects(spxp, phy)
divLeinster(chao.objs$pi, chao.objs$Z, q = 2)
}
