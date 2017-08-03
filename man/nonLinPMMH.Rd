\name{nonLinPMMH}
\alias{nonLinPMMH}
\title{Particle marginal Metropolis-Hastings for a non-linear state space model.}
\usage{
nonLinPMMH(data, particles = 5000, iterations = 10000, burnin = 0,
  plot = FALSE)
}
\arguments{
\item{data}{A vector of the observed data.}

\item{particles}{An integer specifying the number of particles in the particle
filtering estimates of the likelihood.}

\item{iterations}{An integer specifying the number of MCMC iterations.}

\item{burnin}{The number of iterations to remove from the beginning of the MCMC chain
(for plotting purposes only).}

\item{plot}{A boolean variable to determine whether to plot the posterior estimates
and MCMC chain.}
}
\value{
A \code{data.frame} containing the chain
of simulated \eqn{\sigma_v} and \eqn{\sigma_w} values, as well as the
corresponding log likelihood estimates and log prior values.
}
\description{
The \code{nonLinPMMH} function implements particle marginal Metropolis Hastings
for the non-linear state space model described in Section 3.1 of
Andrieu et al. (2010).
}
\details{
This example uses particle marginal Metropolis Hastings to estimate
the standard deviation of the evolution and observation noise in the following
non-linear state space model:

	\eqn{x(n) = 0.5 x(n-1) + 25 x(n-1) / (1+x(n-1)^2) + 8 cos(1.2n)+ e(n)} and 

	\eqn{y(n) = x(n)^2 / 20 + f(n)}

where e(n) and f(n) are mutually-independent normal random
variables of variances var_evol and var_obs, respectively,
and \eqn{x(0) ~ N(0,5)}.

Following Andrieu, Doucet and Holenstein (2010), the priors are
\eqn{var_evol ~ IG(0.01,0.01)} and \eqn{var_obs ~ IG(0.01,0.01)} where IG
is the inverse gamma distribution.

Data can be simulated from the model using \code{\link{simNonlin}}.
}
\examples{
\dontrun{
sim <- simNonlin(len=500,var_init=5,var_evol=10,var_obs=1,cosSeqOffset=0)
res <- nonLinPMMH(sim$data,particles=5000,iterations=50000,burnin=10000,plot=TRUE)
}

}
\references{
C. Andrieu, A. Doucet, and R. Holenstein. Particle Markov chain Monte Carlo methods.
Journal of the Royal Statistical Society: Series B (Statistical Methodology),
72(3):269-342, 2010.
}
\seealso{
\code{\link{simNonlin}} for a function to simulate from the model and
\code{\link{pfNonlinBS}} for a simple bootrap particle filter
applied to a similar non-linear state space model.
}
\author{
Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
}
\concept{
Bayesian PMMH PF
}
\keyword{programming}
