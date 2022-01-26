\name{divLeinster}
\alias{divLeinster}
\title{ Evaluating diversity with a parametrization of the dominance effect}
\usage{
divLeinster(spxp, Z=NULL, q=2, check = TRUE)
}
\description{
divLeinster calculates the diversity of each site of a site by species matrix according to the q parameter according to Leinster & Cobbold 2012.
}
\arguments{
\item{spxp}{Site by species matrix}
\item{Z}{Matrix of similarities among species}
\item{q}{Value of the dominance parameter. Defaults to q = 2.}
\item{check}{Should the arguments be checked?}
}
\references{
Chalmandrier et al. (2015) Effects of species' similarity and dominance on the functional and phylogenetic structure of a plant meta‐community. Ecology 96, 143-153                    Leinster & Cobbold (2012) Measuring diversity: the importance of species similarity. Ecology 93, 477-489
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

divLeinster(spxp, Z, q = 0)
divLeinster(spxp, Z, q = 1)
}
