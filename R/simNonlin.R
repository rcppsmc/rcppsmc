#' Simulates from the associated model.
#' 
#' The \code{simNonlin} function simulates data from the associated model. 
#' 
#' @param len The length of data sequence to simulate.
#' 
#' @return The \code{simNonlin} function returns a list containing the state and data sequences. 
#' 
#' @details The \code{simNonlin} function simulates from the same model.
#' 
#' @rdname pfNonlinBS
#' @export
simNonlin <- function(len = 50)
{
   sim <- list()

   innovations <- rnorm(len) * sqrt(10)
   sim$state[1] <- innovations[1]
   for (i in 2:len) {
       sim$state[i] <- 0.5 * sim$state[i-1] + 25 * sim$state[i-1] /
           (1 + sim$state[i-1]^2) + 8 * cos(1.2*(i-1)) + innovations[i]
   }
   sim$data <- sim$state^2 / 20 + rnorm(len)

   invisible(sim)
}
