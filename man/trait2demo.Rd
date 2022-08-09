\name{Trait2Demo}
\alias{Trait2Demo}
\title{Link function between traits and demography}
\usage{
Trait2Demo(tr, phis, a,b, exp = F)
}
\description{
Allows to calculate a demographic rate based on a linear function of N functional traits. The relationship is determined by N-1 angle coefficients.
}
\arguments{
\item{tr}{Matrix of functional traits. Species as rows, traits as columns.}
\item{phis}{Vectors of phi coefficients}
\item{a}{The variance of the trait}
\item{b}{The average value of the trait}
\item{exp}{is the relationship log-linear ?}
}
\references{
Chalmandrier et al. (2021) Linking functional traits and demography to model species-rich communities. Nature Communications 12, 1-9
}
\examples{
phis <- c(pi/4, 0.56*pi, 4*pi/3)
tr <- matrix(rnorm(40), ncol = 4)

print(Trait2Demo(tr, phis, a = 0.5,b = 1, exp = TRUE))
}
