#' Simulates from the associated model.
#' 
#' The \code{simGaussian} function simulates data from the associated linear
#' Gaussian state space model.
#' 
#' @param len The length of the data sequence to simulate.
#' 
#' @return The \code{simGaussian} function returns a list containing the state and data
#' sequences.
#' 
#' @details The \code{simGaussian} function simulates from the same model returning both
#' the state and observation vectors.
#' 
#' @rdname blockpfGaussianOpt
#' @export
simGaussian <- function(len = 250)
{
  sim <- list()
  sim$state <- cumsum(rnorm(len))
  sim$data  <- sim$state + rnorm(len)

  invisible(sim)
}
