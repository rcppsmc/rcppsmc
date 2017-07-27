\name{blockpfGaussianOpt}
\alias{blockpfGaussianOpt}
\alias{simGaussian}

\title{Block Sampling Particle Filter (Linear Gaussian Model; Optimal Proposal)}
\description{
  The \code{blockpfGaussianOpt} function provides a simple example for
  \pkg{RcppSMC}. It is based on a block sampling particle filter for a linear
  Gaussian model. This is intended only to illustrate the potential of block
  sampling; one would not ordinarily use a particle filter for a model in
  which analytic solutions are available. The 'optimal' block sampler in the
  sense of Doucet, Briers and Senecal (2006) can be implemented in this case.

  The \code{simGaussian} function simulates data from the associated linear
  Gaussian state space model.
}
\usage{
  blockpfGaussianOpt(data, particles=1000, lag=5, plot=FALSE) 
  simGaussian(len)
}
\arguments{
  \item{data}{A vector variable containing the sequence of observations.}
  \item{particles}{An integer specifying the number of particles.}
  \item{lag}{An integer specifying the length of block to use.}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated path along with the uncertainty.}

  \item{len}{The length of the data sequence to simulate.}
}
\value{
  The \code{blockpfGaussianOpt} function returns a matrix containing the final
  sample paths and a vector containing their weights. The logarithm of the
  estimated ratio of normalising constants between the final and initial
  distributions is also returned.

  The \code{simGaussian} function returns a list containing the state and data
  sequences. 
}
\details{
  The \code{blockpfGaussianOpt} function provides a simple example for
  \pkg{RcppSMC}. It is based on a simple linear Gaussian state space model in
  which the state evolution and observation equations are:
  	x(n) = x(n-1) + e(n) and 
	y(n) = x(n) + f(n)
  where e(n) and f(n) are mutually-independent standard normal random
  variables. The 'optimal' block-sampling proposal described by Doucet
  et al (2006) is employed. 

  The \code{simGaussian} function simulates from the same model returning both
  the state and observation vectors.
}
\references{
  A. Doucet, M. Briers, and S. Senecal. Efficient Block Sampling Strategies
  for sequential Monte Carlo methods. Journal of Computational and Graphical
  Statistics, 15(3):693-711, 2006.
}
\examples{
  sim <- simGaussian(len=250)
  res <- blockpfGaussianOpt(sim$data,lag=5,plot=TRUE)
}
\author{Adam M. Johansen and Dirk Eddelbuettel}
\keyword{programming}
