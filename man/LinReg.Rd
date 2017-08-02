\name{LinReg}
\alias{LinReg}
\alias{LinRegLA}
\title{Simple Linear Regression}
\usage{
LinReg(model, particles = 1000, plot = FALSE)

LinRegLA(model, particles = 1000, temperatures = seq(0, 1, 0.05)^5)
}
\arguments{
\item{model}{Choice of regression model (1 for density as the predictor
and 2 for adjusted density as the predictor).}

\item{particles}{An integer specifying the number of particles.}

\item{plot}{A boolean variable to determine whether to plot the posterior estimates.}

\item{temperatures}{In likelihood annealing SMC the targets are defined as \eqn{P(y|\theta)^{\gamma_t}P(\theta)}
where \eqn{0=\gamma_0\le \ldots \le \gamma_T = 1} can be referred to
as the temperatures, \eqn{P(y|\theta)} is the likelihood and \eqn{P(\theta)}
is the prior.}
}
\value{
The \code{LinReg} function returns a \code{list} containing the final particle
approximation to the target (\eqn{\theta} and the corresponding weights) as well as the logarithm
of the estimated model evidence.

The \code{LinRegLA} function returns a \code{list} containing the population
of particles and their associates log likelihoods, log priors and weights at each iteration. The
effective sample size at each of the iterations and several different
estimates of the logarithm of the model evidence are also returned.
}
\description{
A simple example based on estimating the parameters of a linear regression model
using

* Data annealing sequential Monte Carlo (\code{LinReg}).

* Likelihood annealing sequential Monte Carlo (\code{LinRegLA}).
}
\details{
Williams (1959) considers two competing linear regression models
for the maximum compression strength parallel to the grain for radiata
pine. Both models are of the form

\eqn{y_i = \alpha + \beta (x_i - \bar{x}) + \epsilon_i},

where \eqn{\epsilon_i ~ N(0,\sigma^2)} and \eqn{i=1,\ldots,42}.
Here \eqn{y} is the maximum compression strength in pounds per square
inch. The density (in pounds per cubic foot) of the radiata pine
is considered a useful predictor, so model 1 uses density for \eqn{x}.
Model 2 instead considers the density adjusted for resin content, which
is associated with high density but not with strength.

This example is frequently used as a test problem in model choice
(see for example Carlin and Chib (1995) and Friel and Pettitt (2008)).
We use the standard uninformative normal and inverse gamma priors for this example
along with the transformation \eqn{\phi=log(\sigma^2)} so that all parameters
are on the real line and \eqn{\theta = [\alpha,\beta,\phi]}.
The evidence can be computed using numerical estimation
for both of the competing models. The log evidence is -309.9 for model 1 and
-301.4 for model 2.

The \code{LinReg} function implements a data annealing approach to this example.

The \code{LinRegLA} function implements a likelihood annealing approach to this example.
}
\examples{
res <- LinReg(model=1, particles=1000, plot=TRUE)

res <- LinRegLA(model=1, particles=1000)

}
\references{
B. P. Carlin and S. Chib. Bayesian model choice via Markov chain
Monte Carlo. Journal of the Royal Statistical Society: Series B
(Statistical Methodology). 57(3):473-484, 1995.

N. Friel and A. N. Pettitt. Marginal likelihood estimation via power posteriors.
Journal of the Royal Statistical Society: Series B (Statistical Methodology).
70(3):589-607, 2008.

Williams, E. (1959). Regression analysis. Wiley.
}
\author{
Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
}
\concept{
static Bayesian annealing

path
}
\keyword{programming}