simNonlin <- function(len = 50, var_init = 10, var_evol = 10, var_obs = 1, cosSeqOffset = -1)
{
   sim <- list()

   sim$state[1] <- rnorm(1) * sqrt(var_init)
   innovations <- rnorm(len-1) * sqrt(var_evol)
   for (i in 2:len) {
       sim$state[i] <- 0.5 * sim$state[i-1] + 25 * sim$state[i-1] /
           (1 + sim$state[i-1]^2) + 8 * cos(1.2*(i+cosSeqOffset)) + innovations[i-1]
   }
   sim$data <- sim$state^2 / 20 + rnorm(len) * sqrt(var_obs)

   invisible(sim)
}
