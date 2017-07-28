#' A package for general implementation of SMC algorithms.
#'
#' RcppSMC can be used to implement SMC algorithms ranging from
#' simple particle filters to the SMC samplers of Del Moral, Doucet and Jasra (2006) 
#' within a generic framework. For a general introduction to SMC methods, see
#' Doucet, Freitas and Gordon (2001) and Doucet and Johansen (2008).
#'
#' RcppSMC is based around an Rcpp integration of the SMC template class library (Johansen, 2009).
#'
#' At present, three examples are included:
#'
#' \itemize{
#'      \item \code{\link{pfLineartBS}}: A simple 'vehicle tracking' problem of 100 observations solved
#'      using a \dQuote{bootstrap} particle filter with 1000 particles.
#'      \item \code{\link{pfNonlinBS}}: A univariate nonlinear state space model solved using a simple
#'      \dQuote{bootstrap} particle filter. 
#'      \item \code{\link{blockpfGaussianOpt}}: A linear Gaussian model solved using a block sampling
#'      particle filter.
#'}
#'
#' Further integration and extensions are planned. 
#'
#' @references
#' P. Del Moral, A. Doucet and A. Jasra. Sequential Monte Carlo Samplers. Journal of the Royal
#' Statistical Society: Series B (Statistical Methodology), 68(3):411-436, 2006.
#'
#' A. Doucet, N. De Freitas and N.J. Gordon (editors). SMC Methods in Practice. Springer-Verlag, 2001.
#'
#' A. Doucet and A.M. Johansen. A Tutorial on Particle Filtering and Smoothing: 15 Years Later.
#' In the Oxford Handbook of Nonlinear Filtering, Oxford University Press, 2011.
#' 
#' A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
#' Journal of Statistical Software, 30(6):1-41, April
#' 2009. \url{http://www.jstatsoft.org/v30/i06/paper}
#'
#' @seealso The SMCTC paper and code at \url{http://www.jstatsoft.org/v30/i06/paper}.
"_PACKAGE"
#> [1] "_PACKAGE"