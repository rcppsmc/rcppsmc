simGaussian <- function(len = 250)
{
  sim <- c()
  sim$state <- cumsum(rnorm(len))
  sim$data  <- sim$state + rnorm(len)

  invisible(sim)
}
