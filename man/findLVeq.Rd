\name{findLVeq}
\alias{findLVeq}
\title{Algorithmic search of a LV feasible equilibrium}
\usage{
findLVeq(K, A)
}
\description{
This function finds an feasible equilibrium point for a Lotka-Volterra system through an interative algorithm. It aims to solve the equation:

%\loadmathjax

%\mjsdeqn{\mathbf{K} - \mathbf{A} \mathbf{P}^* = 0}
K - AP* = 0
With K the vector of carrying capacities, K the matrix of pairwise interaction coefficients and P* the vector of species abundances at equilibrium to estimate.

The algorithm does the appropriate matrix inversion to solve the equation. If at least one of the abundance value is negative, then the algorithm assumes that the species with the most negative abundance at equilibrium is absent from the community and its abundance is set to 0 and the operation is repeated among the remaining species. This is iterated until all species abundances at equilibrium are positive.

If the A matrix is not inversible, then the Moore-Penrose inverse is instead calculated using \link[corpcor]{pseudoinverse}.
}
\arguments{
\item{K}{Matrix of species traits. Rows represent species and columns represent traits.}
\item{A}{Value of the maximum height of the bell-shaped curve}
}
\references{
Chalmandrier et al. (2022) Predictions of biodiversity are improved by integrating trait-based competition with abiotic filtering. Ecology Letters.
}
\examples{
set.seed(923)
K <- rlnorm(15)
A <- matrix(runif(15*15), ncol = 15)

findLVeq(K, A)
}
