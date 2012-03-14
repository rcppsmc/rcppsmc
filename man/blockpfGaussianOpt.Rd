\name{blockpfGaussianOpt}
\alias{blockpfGaussianOpt}
\title{Block Sampling Particle Filter (Linear Gaussian Model; Optimal Proposal)}
\description{
  The \code{blockpfGaussianOpt} function provides a simple example for
  \pkg{RcppSMC}. It is based on a block sampling particle filter for a linear
  Gaussian model. This is intended only to illustrate the potential of block
  sampling; one would not ordinarily use a particle filter for a model in
  which analytic solutions are available. The 'optimal' block sampler in the
  sense of Doucet, Briers and Senecal (2006) can be implemented in this case.
}
\usage{
  blockpfGaussianOpt(data, particles=1000, lag=5, plot=FALSE) 
}
\arguments{
  \item{data}{A vector variable containing the sequence of observations.}
  \item{particles}{An integer specifying the number of particles.}
  \item{lag}{An integer specifying the length of block to use.}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated path along with the uncertainty.}
}
\value{
  The function returns a matrix containing the final sample paths and a vector
  containing their weights.
}
\details{
  The \code{pfEx} function provides a simple example for
  \pkg{RcppSMC}. It is based on a simple linear Gaussian state space model in
  which the state evolution and observation equations are:
  	x(n) = x(n-1) + e(n) and 
	y(n) = x(n) + f(n)
  where e(n) and f(n) are mutually-independent standard normal random variables. The 'optimal'
  block-sampling proposal described by (Doucet et al., 2006) is employed.
}
\references{
  A. Doucet, M. Briers, and S. Senecal. Efficient Block Sampling Strategies
  for sequential Monte Carlo methods. Journal of Computational and Graphical
  Statistics, 15(3):693-711, 2006.
}
\examples{
\dontrun{
  res <- blockpfGaussianOpt(data,lag=5,plot=TRUE)
}
}
\author{Adam M. Johansen}
\keyword{programming}
