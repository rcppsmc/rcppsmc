simGaussianSSM <- function(len = 100,
                           stateInit = 0,
                           phi = 0.7,
                           varStateEvol = 1,
                           varObs = 1) {

    sim <- list("states" = rep(0, times = len),
                "measurements"  = rep(0, times = len))
    stateEvol   <- rnorm(len) * sqrt(varStateEvol)
    observErr   <- rnorm(len) * sqrt(varObs)

    if (is.null(stateInit)) {
        stateInit <- rnorm(1, 0, sqrt(varStateEvol / (1 - phi^2)))
    }
    sim$states[1] <- phi * stateInit + stateEvol[1]
    for (i in 2:len) {
        sim$states[i] <- phi * sim$states[i - 1] + stateEvol[i]
    }
    sim$measurements <- sim$states + observErr
    return(sim)
}
