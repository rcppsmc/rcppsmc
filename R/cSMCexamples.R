compareNCestimates <- function(data,
                               numParticles = 1000L,
                               simNumber = 10L,
                               parameters = list(phi = 0.9,
                                                 varStateInit = 10,
                                                 varStateEvol = 10,
                                                 varObs = 1),
                               plot = FALSE) {

    # if no data supplied, use default
    if (missing(data)) data <- simGaussianSSM()

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) != 2)

    browser()
    res <- compareNCestimates_imp(as.vector(data),
                                  numParticles,
                                  simNum = simNumber,
                                  parInits = parameters)

    invisible(res)
}
