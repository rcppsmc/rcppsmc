\name{pfLineartBS}
\alias{pfLineartBS}
\alias{pfExOnlinePlot}
\title{Particle Filter Example}
\description{
  The \code{pfLineartBS} function provides a simple example for
  \pkg{RcppSMC}. It is based on the first example in \code{SMCTC} and
  the discussion in Section 5.1 of Johansen (2009). A simple 'vehicle
  tracking' problem of 100 observations is solved with 1000 particles.

  The \code{pfExOnlinePlot} function provides a simple default
  \sQuote{online} plotting function that is invoked during the
  estimation process. 
}
\usage{
  pfLineartBS(data, particles=1000, plot=FALSE, onlinePlot) 
  pfExOnlinePlot(xm, ym)
}
\arguments{
  \item{data}{A two-column matrix or dataframe containing x and y
    values. The default data set from Johansen (2009) is used as the
    default if no data is supplied.}
  \item{particles}{An integer specifying the number of particles.}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated path along with the data.}
  \item{onlinePlot}{A user-supplied callback function which is called with the
  x and y position vectors during each iteration of the algorithm; see
  pfExOnlinePlot for a simple example.}
  \item{xm}{Vector with x position.}
  \item{ym}{Vector with y position.}
}
\value{
  The function returns a \code{data.frame} containing as many rows as in
  the input data, and four columns corresponding to the estimated \eqn{x}{x} and
  \eqn{y}{y} coordinates as well as the estimated velocity in these two directions.
}
\details{
  The \code{pfLineartBS} function provides a simple example for
  \pkg{RcppSMC}. The model is linear with t-distributed innovations.
  It is based on the \code{pf} example in the
  \code{SMCTC} library, and discussed in the Section 5.1 of his
  corresponding paper (Johansen, 2009).

  Using the simple \code{pfExOnlinePlot} function illustrates how
  callbacks into R, for example for plotting,  can be made during the
  operation of SMC algorithm.
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{http://www.jstatsoft.org/v30/i06/paper}
}
\seealso{The SMCTC paper and code at \url{http://www.jstatsoft.org/v30/i06/paper}.}
\examples{
\dontrun{
  res <- pfLineartBS(plot=TRUE)
  res <- pfLineartBS(onlinePlot=pfExOnlinePlot)
}
}
\author{Adam M. Johansen for SMCTC, Dirk Eddelbuettel for the RcppSMC wrapper.}
\keyword{programming}
