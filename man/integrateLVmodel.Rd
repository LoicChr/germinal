\name{integrateLVmodel}
\alias{integrateLVmodel}
\title{Run a Lotka-Volterra model}
\usage{
integrateLVmodel(r,a,init,t0, t1, dt)
}

\arguments{
\item{r}{Vector of intrinsic growth rates}
\item{a}{Matrix of pairwise interactions coefficients}
\item{init}{Vector of species initial abundances}
\item{t0}{Initial time step}
\item{t1}{Final time step}
\item{dt}{Time step length}
}
\description{
Integrates an ordinary differential equation (ODE) model of a Lotka-Volterra system of species characterized by their intrinsic growth rate and their pairwise interaction coefficients. The ODE follows the following model with :

\loadmathjax

%\mjsdeqn{\frac{dP_i}{dt} = P_i (r_i - \sum_j a_{ij} P_j )}
dPi/dt = Pi (ri - sum aij Pj)

With Pi the abundance of species i, r_i its intrinsic growth rate and a_ij the (negative) per capita impact of species j on species i.


The model uses the c++ library odeint.
}

\references{
Chalmandrier L. et al. (2021) Linking functional traits and demography to model species-rich communities. Nature Communications 12: 2724

Lotka, A. J. (1926) Elements of physical biology. Williams & Wilkins Co, Baltimore
}
\examples{
out = integrateLVmodel(r = c(1,2),a = matrix(c(0.5,-0.1,-1.2,0.8), ncol = 2),c(7,7),t0 = 0, 3, 0.01)
plot(out[,1], out[,2], ylim = c(0,15), type ="b", xlab = "Time", ylab = "Species abundance")
points(out[,1], out[,3], type = "b", col = "red")
}

