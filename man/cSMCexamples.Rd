\name{compareNCestimates}
\alias{compareNCestimates}
\alias{simGaussianSSM}
\alias{kalmanFFBS}
\title{Conditional Sequential Monte Carlo Examples}
\description{
  The \code{compareNCestimates} function generates a Monte Carlo study to
  compare log-likelihood (normalizing constant) estimates in the standard linear
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
                                            phi = 0.7,
                                            varStateEvol = 1,
                                            varObs = 1),
                     plot = FALSE)
  simGaussianSSM(len = 100,
                 stateInit = 0,
                 phi = 0.7,
                 varStateEvol = 1,
                 varObs = 1)
  kalmanFFBS(dataY,
             stateInit,
             phi,
             varStateEvol,
             varObs,
             simNumber)
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
  \code{compareNCestimates} returns a named list of two
  \itemize{
    \item{\code{outSMC}} a named list of two:
    \itemize{
      \item{\code{smcOut}: }{a matrix of dimension \code{simNum x 4} which
      contains single log-likelihood estimates of the standard BPF for each of
      the 4 resampling schemes and for each simulation run}
      \item{\code{csmcOut}: }{a matrix of dimension \code{simNum x 4} which
      contains single log-likelihood estimates of the conditional BPF for each
      of the 4 resampling schemes and for each simulation run}
    }
    \item{\code{outKalman}} the output of \code{kalmanFFBS}, see below
  }
  \code{kalmanFFBS} returns a named list of two:
    \itemize{
      \item{\code{logLikeliEstim: }}{(exact) estimate of the log-likelihood}
      \item{\code{xBackwardSimul: }}{a backward simulation (joint smoothing) trajectory of latent states given parameters and measurement}
    }
}
\details{
  \code{compareNCestimates} runs a simulation study that provides
  log-likelihood (normalizing constant) estimates; there are \code{simNumber}
  runs of the standard BPF and the conditional BPF under four resampling schmes:
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

  \code{simGaussianSSM} simulates from a Linear Gaussian state space model of
  the following form:
  \deqn{x_t=\phi x_{t-1} + u_t}
  \deqn{y_t= x_{t} + w_t}
  where \eqn{\phi} is set via the \code{phi} argument,
  \eqn{u_t \sim N(0, \sigma_x^2)}, \eqn{w_t \sim N(0, \sigma_y^2)} for which the
  innovation (\eqn{\sigma_x^2}) and measurement (\eqn{\sigma_y^2}) variances are
  set via arguments \code{varStateEvol} and \code{varObs}, respectively.
}
\references{
  A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
  Journal of Statistical Software, 30(6):1-41, April
  2009. \url{https://www.jstatsoft.org/article/view/v030i06}
}

\seealso{The SMCTC paper and code at \url{https://www.jstatsoft.org/article/view/v030i06}.}

\author{Adam M. Johansen, Dirk Eddelbuettel, Leah South and Ilya Zarubin}
\keyword{programming}
