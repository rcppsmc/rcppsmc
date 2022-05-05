\name{compareNCestimates}
\alias{compareNCestimates}
\alias{simGaussianSSM}
\alias{kalmanFFBS}
\title{Conditional Sequential Monte Carlo Examples}
\description{
  The \code{compareNCestimates} function generates a Monte Carlo study to
  compare likelihood (normalizing constant) estimates in the standard linear
  Gaussian state space (LGSS) model: Kalman filter estimates, as the benchmark,
  are compared to the standard bootstrap particle filter and the conditional
  bootstrap particle filter estimates (see \code{Details}).


  The \code{simGaussianSSM} function simulates data from a LGSS model (can be
  used manually to simulate \code{data} or runs as a default, if no \code{data} is provided, with a default parameter setup, see \code{parameters}).

  The \code{kalmanFFBS} function runs a Kalman (exact) forward filter, computes
  a log-likelihood estimate and generates a joint smoothing state trajectory
  via a backward simulation pass.
}
\usage{
  compareNCestimates(dataY,
                     trueStates = NULL,
                     numParticles = 1000L,
                     simNumber = 100L,
                     modelParameters = list(stateInit = 0,
                                            phi = 0.9,
                                            varStateEvol = 1,
                                            varObs = 1),
                     plot = FALSE)
  simGaussianSSM(len = 100,
                 stateInit = 0,
                 phi = 0.9,
                 varStateEvol = 1,
                 varObs = 1)
  kalmanFFBS(dataY,
             stateInit,
             phi,
             varStateEvol,
             varObs)
}
\arguments{
  \item{dataY}{A one-column matrix or dataframe or vector containing
  measurements (y values) from a standard linear Gaussian SSM. If not provided, defaults to a LGSS model with time series lenght \code{len=250} and parameter
  setup specified with default values in the \code{parameters} argument, see \code{\link{simGaussianSSM}} or \code{\link{compareNCestimates}}.}
  \item{trueStates}{defaults to \code{NULL} for a real dataset as the true state
  values are not observed; for simulated data, these can be passed and then will alse be plotted if \code{plot=TRUE}}
  \item{numParticles}{An integer specifying the number of particles.}
  \item{simNumber}{An integer specifying the number of repeated simulation runs of each of which produces 2x4=8 normalizing constant esimtates: BPF and conditional BPF esimates under four conditional resampling schemes, as well as
  a ground truth Kalman forward filter estimate and a backward filter output
  required for the reference trajectory of the conditional sampler.}
  \item{modelParameters}{a named \code{list} of parameters of the LGSSM model
  in the following order:
  \itemize{
    \item \code{phi:} autoregressive parameter
    \item \code{stateInit:} initial state value (i.e. \eqn{X_0})
    \item \code{varStateEvol:} state process variance
    \item \code{varObs:} measurement/observation process variance
  }
  These parameters are used to for the Kalman forward filtering and backward
  simulation pass, and, if no \code{data} argument is provided, to simulate
  data from a linear Gaussian state space model internally via
  \code{simGaussianSSM}.}
  \item{phi}{autoregressive parameter}
  \item{stateInit}{initial state value (i.e. \eqn{X_0})}
  \item{varStateEvol}{state process variance}
  \item{varObs}{measurement/observation process variance}
  \item{plot}{A boolean variable describing whether plot should
    illustrate the estimated results along with the data.}
  \item{len}{Length of data series to simulate.}
}
\value{
  \code{compareNCestimates} function returns a \code{data.frame} containing as many rows as in
  the input data, and four columns corresponding to the estimated \eqn{x}{x} and
  \eqn{y}{y} coordinates as well as the estimated velocity in these two
  directions.

  \code{simLineart} function returns a list containing the vector of
  states and the associated vector of observations.

  \code{kalmanFFBS} returns a named list of two:
    \itemize{
      \item{\code{logLikeliEstim:}}{(exact) estimate of the log-likelihood}
      \item{\code{xBackwardSimul:}}{a backward simulation (joint smoothing) trajectory of latent states given parameters and measurement}
    }
}
\details{
  The \code{compareNCestimates} The estimates are obtained from \code{sumNumber} runs of the standard BPF and the conditional BPF under four
  resampling schmes:
  \itemize{
    \item multinomial
    \item stratified
    \item systematic
    \item residual
  }
  The "ground truth" Kalman forward filter estimate of the normalizing constant
  is compared to the BPF normalizing constant estimates, which are unbiased for
  all above schemes; specifically the conditional BPF estimate is unbiased if
  the reference trajectory is simulated from the target distribution which is obtained here as a backward simulation run of the Kalman
  filter.

  Box plots illustrate the unbiasedness of standard BPF and conditional BPF
  estimates for the normailizing constant estimate in the linear Gaussian SSM, and serve as an small example for to illustrate conditional SMC
  algorithms (in their most basic BPF version) with different conditonal resampling schemes as implemented within \pkg{RcppSMC}.
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{https://www.jstatsoft.org/v30/i06/paper}
}

\seealso{The SMCTC paper and code at \url{https://www.jstatsoft.org/v30/i06/paper}.}

\author{Adam M. Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin}
\keyword{programming}
