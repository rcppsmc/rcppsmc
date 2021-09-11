\name{compareNCestimates}
\alias{compareNCestimates}
\alias{simGaussianSSM}
\title{Conditional Sequential Monte Carlo Examples}
\description{
  The \code{compareNCestimates} The function generates a Monte Carlo study to
  compare likelihood (normalizing constant) estimates for the standard linear
  Gaussian state space model.

  The \code{simGaussianSSM} function simulates data from the linear Gaussian state space model.
}
\usage{
  compareNCestimates(data,
                     particles = 1000L,
                     simNumber = 10L,
                     parameters = list(phi = 0.9,
                                       varInit = 10,
                                       varEvol = 10,
                                       varObs = 1),
                     plot = FALSE)
  simGaussianSSM(len = 250,
                 parameters = list(phi = 0.9,
                                   varInit = 10,
                                   varEvol = 10,
                                   varObs = 1))
}
\arguments{
  \item{data}{A one-column matrix or dataframe or vector containing measurements
   (y values). Defaults to a standard linear Gaussian SSM with time series lenght 250 (for defaults paraemter settings, see function help.)}
  \item{particles}{An integer specifying the number of particles.}
  \item{simNumber}{An integer specifying the number of repeated simulation runs of each of which produces 6 normalizing constant esimtates: BPF, four conditional SMC estimates corresponding to four different resampling schemes,
  and a ground truth Kalman filter estimate.}
  \item{parameters}{a \code{list} containing named parameters of the linear
  state space model in the following order: phi (autoregressive parameter), varInit (initial variance), varEvol (state process variance), and varObs (
      measuremnt variance
  )}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated path along with the data.}
  \item{len}{Length of data sequence to simulate.}
}
\value{
  The \code{compareNCestimates} function returns a \code{data.frame} containing as many rows as in
  the input data, and four columns corresponding to the estimated \eqn{x}{x} and
  \eqn{y}{y} coordinates as well as the estimated velocity in these two
  directions.

  The \code{simLineart} function returns a list containing the vector of
  states and the associated vector of observations.
}
\details{
  The \code{compareNCestimates} The estimates are taken from the standard BPF,
  a conditional SMC with different resampling schemes: multinomial, stratified,
  systematic and residuals. The normalizing constant estimates are compared to
  the "ground trouth" Kalman Filter (backward smoothing) output.

  This function also provides a small example for implementing conditional SMC
  algorithms with different conditonal resampling schemes implemented currently
  in \pkg{RcppSMC}.
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{https://www.jstatsoft.org/v30/i06/paper}
}
# \examples{
#   data <- simGaussianSSM(len = 250,
#                          parameters = list(phi = 0.9,
#                                            varInit = 10,
#                                            varEvol = 10,
#                                            varObs = 1))
#   res  <- compareNCestimates(data,
#                                particles = 1000L,
#                                simNumber = 10L,
#                                parameters = list(phi = 0.9,
#                                                  varInit = 10,
#                                                  varEvol = 10,
#                                                  varObs = 1),
#                                plot = FALSE)
# }
\seealso{The SMCTC paper and code at \url{https://www.jstatsoft.org/v30/i06/paper}.}

\author{Adam M. Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin}
\keyword{programming}
