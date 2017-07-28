#' Simulates from the associated model.
#' 
#' The \code{simLineart} function simulates data from the model. 
#' 
#' @param len Length of sequence to simulate.
#' 
#' @return The \code{simLineart} function returns a list containing the vector of
#' states and the associated vector of observations.
#' 
#' @details \code{simLineart} simulates from the model. 
#' 
#' @rdname pfLineartBS
#' @export
simLineart <- function(len = 250)
{
  sim <- list()

  statex <- cumsum(rnorm(len))
  statey <- cumsum(rnorm(len))

  datax  <- statex + rt(len,df=4)
  datay  <- statey + rt(len,df=4)

  sim$state <- matrix(cbind(statex,statey), ncol=2, dimnames=list(NULL, c("x", "y")))
  sim$data  <- matrix(cbind(datax,datay), ncol=2, dimnames=list(NULL, c("x", "y")))

  invisible(sim)
}
