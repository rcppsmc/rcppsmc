\name{pfNonlinBS}
\alias{pfNonlinBS}
\title{Block Sampling Particle Filter (Linear Gaussian Model; Optimal Proposal)}
\description{
  The \code{blockpfGaussianOpt} function provides a simple example for
  \pkg{RcppSMC}. It is a simple "bootstrap" particle filter which employs
  multinomial resampling after each iteration applied to the common nonlinear
  state space model following Gordon, Salmond and Smith (1993).
}
\usage{
  pfNonlinBS(data, particles=500, plot=FALSE) 
}
\arguments{
  \item{data}{A vector variable containing the sequence of observations.}
  \item{particles}{An integer specifying the number of particles.}
  \item{plot}{A boolean variable describing whether a plot should
    illustrate the (posterior mean) estimated path along with one and two
    standard deviation intervals.} 
}
\value{
  The function returns two vectors, the first containing the posterior
  filtering means; the second the posterior filtering standard deviations.
}
\details{
  The \code{pfNonlinbs} function provides a simple example for
  \pkg{RcppSMC}. It is based on a simple nonlinear state space model in
  which the state evolution and observation equations are:
  	x(n) = 0.5 x(n-1) + 25 x(n-1) / (1+x(n-1)^2) + 8 cos(1.2(n-1))+ e(n) and 
	y(n) = x(n)^2 / 20 + f(n)
  where e(n) and f(n) are mutually-independent normal random
  variables of variances 10.0 and 1.0, respectively. A boostrap proposal
  (i.e. sampling from the state equation) is used, together with multinomial
  resampling after each iteration. 
}
\references{
  N. J. Gordon, S. J. Salmond, and A. F. M. Smith. Novel approach to
  nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings-F, 
  140(2):107-113, April 1993.	
}
\examples{
\dontrun{
  res <- pfNonlinBS(data,particles=500,plot=TRUE)
}
}
\author{Adam M. Johansen}
\keyword{programming}
