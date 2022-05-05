\name{RcppSMC-Functions}
\alias{rcppsmc_logNormWeightsCpp}
\title{Example functions for computing stable log sum weights}
\description{
  This function is created when
  \code{RcppSMC.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcppsmc_logNormWeightsCpp(logWeights)
}
\arguments{
  \item{logWeights}{a numeric vector with logarithmic weights}
}
\value{
  \code{rcppsmc_logNormWeightsCpp()} returns a named list of two, the
  logarithmic sum of weights and the normalised logarithmic weights.
}
\details{
  This is an example function which computes stable calculations of the
  logarithm of the sum of weights (which is of interest for computing e.g.
  normalising constants) as well as the normalised logarithmic weights. This is
  to demonstrate how to write a function using the RcppSMC C++ internal
  functions and classes, and to have to such a function accessible from R.
}
\references{
  See the documentation for RcppSMC for more details.
}
\examples{
  logWeightVec <- -1 * sqrt(1:4)
  rcppsmc_logNormWeightsCpp(logWeightVec)
}
\author{Dirk Eddelbuettel, Adam Johansen, Leah F. South, Ilya Zarubin}

