\name{rareEventsEx}
\alias{rareEventsEx}
\title{Particle Filter Example}
\description{
  The \code{rareEventsEx} function provides another example for
  \pkg{RcppSMC}. It is based on the second example in \code{SMCTC} and
  the discussion in Section 5.2 of Johansen (2009). It implements an
  algorithm in Johansen et al (2006) which estimates rare event
  probabilities for a Gaussian distribution.
}
\usage{
  rareEventsEx(number=100, iterations=10, threshold=5.0, schedule=30.0)
}
\arguments{
  \item{number}{An integer specifying the number of particles.}
  \item{iterations}{An integer specifying the number of iterations.}
  \item{threshold}{A numericd describing the rare event threshold.}
  \item{schedule}{A numericd describing the annealing schedule constant.}
}
\value{
  The function returns a named \code{vector.frame} containing three values.
}
\details{
  The \code{rareEventsEx} function provides the second simple example for
  \pkg{RcppSMC}. It is based on the \code{rare-events} example in the
  \code{SMCTC} library, and discussed in the Section 5.2 of his
  corresponding paper (Johansen, 2009), which in turn refers to Johansen
  et al (2006) for a deeper discussion.
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{http://www.jstatsoft.org/v30/i06/paper}

  A. M. Johansen, P. Del Moral, A. Doucet. Sequential Monte Carlo
  Samplers for Rare Events. Proceedings of the 6th International
  Workshop on Rare Event Simulation, pages 256--267, Bamberg (Germany).
}
\seealso{The SMCTC paper and code at \url{http://www.jstatsoft.org/v30/i06/paper}.}
\examples{
\dontrun{
  res <- rareEventsEx()
}
}
\author{Adam M. Johansen for SMCTC, Dirk Eddelbuettel for the RcppSMC wrapper.}
\keyword{programming}
