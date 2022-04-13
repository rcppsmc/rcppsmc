simGaussianSSM <- function(len = 250,
                           parameters = list(phi = 0.9,
                                             varInit = 1,
                                             varEvol = 1,
                                             varObs = 1)) {


    sim <- list()
    sim$state[1] <- rnorm(1) * sqrt(parameters$varInit)
    innovations  <- rnorm(len - 1) * sqrt(parameters$varEvol)
    for (i in 2:len) {
        sim$state[i] <- parameters$phi * sim$state[i - 1] + innovations[i - 1]
    }
    sim$data <- sim$state + rnorm(len) * sqrt(parameters$varObs)
    return(sim)
}
