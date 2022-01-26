\name{abgDecompQ}
\alias{abgDecompQ}
\title{ Alpha, beta, gamma decomposition with parametrization of the dominance effect}
\usage{
abgDecompQ(spxp, Z=NULL, q=2, check=TRUE)
}
\description{
divLeinster calculates the diversity of each site of a site by species matrix according to the q parameter according to Leinster & Cobbold 2012.
abgDecompQ performs a alpha, beta, gamma multiplicative decomposition using Leinster's diversity indices.

BetaDisQ calculates the pairwise beta-diversity (minus 1) between sites of a site by species matrix according to the q parameter using the afformentionned functions
Allows a parametrization of the dominance effect

chaoObjects is a data preparation function. It returns adequate arguments for abgDecompQ, BetaDisQ and divLeinster to perform a diversity analysis using Chao's diversity index.
Warning: this formula only works with an ultrametric tree!
}
\arguments{
\item{spxp}{Site by species matrix}
\item{Z}{Matrix of similarities among species}
\item{q}{Value of the dominance parameter. Defaults to q = 2.}
\item{check}{Should the arguments be checked?}
}
\references{
Chalmandrier et al. (2015) Effects of species' similarity and dominance on the functional and phylogenetic structure of a plant meta‐community. Ecology 96, 143-153                    Chao et al. 2010 Phylogenetic diversity measures based on Hill numbers. Philosophical Transactions of the Royal Society B: Biological Sciences 365: 3599-3609                          Leinster & Cobbold (2012) Measuring diversity: the importance of species similarity. Ecology 93, 477-489
}
\examples{
library(ape)
data(phylocom, package = "picante")

spxp <- phylocom$sample
phy <- phylocom$phylo
is.missing<- !(phy$tip.label \%in\% colnames(spxp))
phy <- drop.tip(phy, phy$tip.label[is.missing])
spxp <- spxp[, phy$tip.label]

dis <- cophenetic(phy)
Z <- 1 - dis/max(dis)

abgDecompQ(spxp, Z, q = 0)
abgDecompQ(spxp, Z, q = 1)
}
