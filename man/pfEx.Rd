\name{pfEx}
\alias{pfEx}
\title{Particle Filter Example}
\description{
  The \code{pfEx} function provides a simple example for
  \pkg{RcppSMC}. It is based on the first example in \code{SMCTC} and
  the discussion in Section 5.1 of Johansen (2009). A simple 'vehicle
  tracking' problem of 100 observations is solved with 1000 particles.
}
\usage{
  pfEx(filename="", particles=1000, plot=FALSE) 
}
\arguments{
  \item{filename}{A character variable describing path and filename of a
    data file in the format described in the \pkg{SMCTC}
    documentation. Note we may replace this with a direct transfer of a
    data matrix.}
  \item{particles}{An integer specifying the number of particles.}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated path along with the data.}
}
\value{
  The function returns a \code{data.frame} containing as many rows as in
  the input data, and four columns corresponding to the estimated \eqn{x}{x} and
  \eqn{y}{y} coordinates as well as the estimated velocity in these two directions.
}
\details{
  The \code{pfEx} function provides a simple example for
  \pkg{RcppSMC}. It is based on the \code{pf} example in the
  \code{SMCTC} library, and discussed in the Section 5.1 of his
  corresponding paper (Johansen, 2009).
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{http://www.jstatsoft.org/v30/i06/paper}
}
\seealso{The SMCTC paper and code at \url{http://www.jstatsoft.org/v30/i06/paper}.}
\examples{
\dontrun{
  res <- pfEx(plot=TRUE)
}
}
\author{Adam M. Johansen for SMCTC, Dirk Eddelbuettel for the RcppSMC wrapper.}
\keyword{programming}
