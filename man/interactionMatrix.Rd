\name{interactionMatrix}
\alias{interactionMatrix}
\title{Trait-dependent pairwise interaction matrix}
\usage{
interactionMatrix(tr, intercept, mu, sigma, rho = NULL, intra = 1, std = F)
}
\description{
Generates a pairwise interaction matrix that assumes that coefficients are a bell-shaped (Gaussian) function of species trait differences.
}
\arguments{
\item{tr}{Matrix of species traits. Rows represent species and columns represent traits.}
\item{intercept}{Value of the maximum height of the bell-shaped curve}
\item{mu}{vector that contains the peak position along the trait difference dimensions.}
\item{sigma}{vector that contains the width of the Gaussian function along the trait difference dimensions}
\item{rho}{vector of partial correlation parameters that specify the interaction among trait difference dimensions. Each element must be between -1 and 1}
\item{intra}{Value of the intraspecific interaction parameters. If NULL, the value is calculated similarly to interspecific interaction coefficients}
\item{std}{Do we standardize the interaction coefficients so that aij(dt = 0) = C}

}
\references{
Chalmandrier et al. (2022) Predictions of biodiversity are improved by integrating trait-based competition with abiotic filtering. Ecology Letters.
}
\examples{
tr <- matrix(rnorm(3*15), ncol = 3)

# Single trait interaction matrix
interactionMatrix(tr = tr[,1, drop = FALSE], intercept = 0.5, mu = -1, sigma =2,  intra =1)

# Three trait interaction matrix
interactionMatrix(tr = tr, intercept = 0.5, mu = c(-1,0,0.5),
                  sigma =c(2,0.5, 6), rho = c(0,0.5,-0.6), intra =1)

}
