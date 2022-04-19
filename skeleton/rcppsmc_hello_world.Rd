\name{RcppArmadillo-Functions}
\alias{rcppsmc_logNormmWeightsCpp}
\title{Example functions for computing stable log sum weights}
\description{
  This function is created when
  \code{RcppSMC.package.skeleton()} is invoked to create a
  skeleton packages.
}
\usage{
rcppsmc_logNormmWeightsCpp()
}
\arguments{
  \item{logWeights}{a numeric vector with logarithmic weights}
}
\value{
  \code{rcppsmc_logNormmWeightsCpp()} returns a named list of two, the
  logarithmic sum of weights and the logarithmic weights as passed via the
  function argument \code{logWeights}.
}
\details{
  This is an example function which computes stable calculations of the
  logarithm of the sum of weights (which is of interest for computing e.g.
  normalising constants). This is to demonstrate how to write a function using
  the RcppSMC C++ internal functions and classes, and to have to such a function
  accessible from R.
}
\references{
  See the documentation for RcppSMC for more details.
}
\examples{
  x <- sqrt(1:4)
  rcpparma_innerproduct(x)
  rcpparma_outerproduct(x)
}
\author{Dirk Eddelbuettel, Adam Johansen, Leah F. South, Ilya Zarubin}

