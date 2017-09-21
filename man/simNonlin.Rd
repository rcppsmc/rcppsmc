\name{simNonlin}
\alias{simNonlin}
\title{Simulates from a simple nonlinear state space model.}
\usage{
simNonlin(len = 50, var_init = 10, var_evol = 10, var_obs = 1,
  cosSeqOffset = -1)
}
\arguments{
\item{len}{The length of data sequence to simulate.}

\item{var_init}{The variance of the noise for the initial state.}

\item{var_evol}{The variance of the noise for the state evolution .}

\item{var_obs}{The variance of the observation noise.}

\item{cosSeqOffset}{This is related to the indexing in the
cosine function in the evoluation equation. A value of -1
can be used to follow the specification of Gordon, Salmond
and Smith (1993) and 0 can be used to follow
Andrieu, Doucet and Holenstein (2010).}
}
\value{
The \code{simNonlin} function returns a list containing the state and data sequences.
}
\description{
The \code{simNonlin} function simulates data from the models used
in \code{link{pfNonlinBS}} and \code{link{nonLinPMMH}}.
}
\details{
The \code{simNonlin} function simulates from
a simple nonlinear state space model with
state evolution and observation equations:
	
	\eqn{x(n) = 0.5 x(n-1) + 25 x(n-1) / (1+x(n-1)^2) + 8 cos(1.2(n+cosSeqOffset))+ e(n)} and 
	
	\eqn{y(n) = x(n)^2 / 20 + f(n)}
	
where \eqn{e(n)} and \eqn{f(n)} are mutually-independent normal random
variables of variances var_evol and var_obs, respectively,
and \eqn{x(0) ~ N(0,var_init)}.

Different variations of this model can be found in
Gordon, Salmond and Smith (1993) and
Andrieu, Doucet and Holenstein (2010). A cosSeqOffset
of -1 is consistent with the former and
0 is consistent with the latter.
}
\references{
C. Andrieu, A. Doucet, and R. Holenstein. Particle Markov chain Monte Carlo methods.
Journal of the Royal Statistical Society: Series B (Statistical Methodology),
72(3):269-342, 2010.

N. J. Gordon, S. J. Salmond, and A. F. M. Smith. Novel approach to
nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings-F, 
140(2):107-113, April 1993.
}
\author{
Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
}
\seealso{
\code{\link{pfNonlinBS}} for a simple bootrap particle filter
applied to this model and \code{\link{nonLinPMMH}} for particle
marginal Metropolis Hastings applied to estimating the standard
deviation of the state evolution and observation noise.
}
