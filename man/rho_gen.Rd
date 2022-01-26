\name{rho_gen}
\alias{rho_gen}
\title{Parametrized correlation matrix}
\usage{
rho_gen(thetas)
}
\description{
\loadmathjax
Generates a correlation matrix among N vectors using \mjseqn{\frac{N(N-1)}{2}} \mjseqn{\theta} parameters. The cosinus of \mjseqn{\theta} parameters are the partial correlation coefficients among vectors.
}
\arguments{
\item{thetas}{vector of \mjseqn{\theta} parameters. Each element must be between 0 and \mjseqn{\pi}}

}
\references{
Forrester, P.J & Zhang, J. (2020) Parametrising correlation matrices. Journal of Multivariate Analysis 178, 104619
}
\examples{
thetas <- runif(6, 0,pi)

rho_gen(thetas)
}
