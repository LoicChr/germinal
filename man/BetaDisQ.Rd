\name{BetaDisQ}
\alias{BetaDisQ}
\title{Pairwise beta-diversity with a parametrization of the dominance effect}
\usage{
BetaDisQ(spxp, Z=NULL, q=2, check = TRUE)
}
\description{
BetaDisQ calculates the pairwise beta-diversity (minus 1) between sites of a site by species matrix according to the q parameter using divLeinster and abgDecompQ
}
\arguments{
\item{spxp}{Site by species matrix}
\item{Z}{Matrix of similarities among species}
\item{q}{Value of the dominance parameter. Defaults to q = 2.}
\item{check}{Should the arguments be checked?}
}
\references{
Chalmandrier et al. (2015) Effects of species' similarity and dominance on the functional and phylogenetic structure of a plant meta‐community. Ecology 96, 143-153
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

BetaDisQ(spxp, Z, q = 0)
BetaDisQ(spxp, Z, q = 1)
}
