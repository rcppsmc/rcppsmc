compareNCestimates <- function(dataY,
                               trueStates = NULL,
                               numParticles = 1000L,
                               simNumber = 100L,
                               modelParameters = list(stateInit = 0,
                                                      phi = 0.9,
                                                      varStateEvol = 1,
                                                      varObs = 1),
                               plot = FALSE) {
    # more eloquent tests can be added
    # stopifnot(nrow(data) > 0,
    #           ncol(data) != 2)

    # if no data supplied, use default
    if (missing(dataY)) {
        dataSim  <- simGaussianSSM()
        dataY    <- dataSim$measurements
        # for simulated data, the true states are available
        trueStates <- dataSim$states
    }

    browser()
    resFFBS <- kalmanFFBS(as.vector(dataY),
                          stateInit = modelParameters$stateInit,
                          phi = modelParameters$phi,
                          varStateEvol = modelParameters$varStateEvol,
                          varObs = modelParameters$varObs)
    resSMC <- compareNCestimates_imp(as.vector(dataY),
                                     numParticles,
                                     simNum = simNumber,
                                     parInits = modelParameters,
                                     resFFBS$xBackwardSimul)
    if (isTRUE(plot)) {
      par(mfrow = c(3, 2))
      if (!is.null(trueStates)) {
          layoutMatrix <- matrix(c(1, 1, 2, 2,
                                   3, 4, 5, 6,
                                   7, 0, 0, 0), nrow = 3, ncol = 4,
                                   byrow = TRUE)

          layout(mat = layoutMatrix,
                 heights = c(1, 1, 1), # Heights of the two rows
                 widths = c(1, 1, 1, 1)) # Widths of the two columns

        title1 <- "measurements (green) vs. true latent states (red)"
        title2 <- "true latent states (red) vs. backward simulation path (blue)"
        plot(data$measurements, type = "l", col = "blue", main = title1)
        lines(trueStates, type = "l", col = "red")
        plot(trueStates, type = "l", col = "red", main = title2)
        lines(resFFBS$xBackwardSimul, type = "l", col = "blue")
      } else {
        title3 <- "measurements (green) vs. backward simulation path (blue)"
        plot(data$measurements, type = "l", col = "blue", main = title3)
        lines(resFFBS$xBackwardSimul, type = "l", col = "red")
        plot.new()
      }
      boxplot(resSMC$SMC$smcOut[, 1])
      abline(h = resFFBS$Kalman$logLikeliEstim, col = "red")
      boxplot(resSMC$SMC$smcOut[, 2])
      abline(h = resFFBS$Kalman$logLikeliEstim, col = "red")
      boxplot(resSMC$SMC$smcOut[, 3])
      abline(h = resFFBS$Kalman$logLikeliEstim, col = "red")
      boxplot(resSMC$SMC$smcOut[, 4])
      abline(h = resFFBS$Kalman$logLikeliEstim, col = "red")
      boxplot(resSMC$SMC$csmcOut[, 1])
      abline(h = resFFBS$Kalman$logLikeliEstim, col = "red")
    }


    # invisible(res)
    return(list(SMC = resSMC, Kalman = resFFBS))
}
kalmanFFBS <- function(data,
                       stateInit,
                       phi,
                       varStateEvol,
                       varObs) {
    y   <- data
    len <- length(y)
    # KF part:
    # Housekeeping:
    ## Forward Filtering Container:
    xt  <- numeric(len)
    Pt  <- numeric(len)
    xt1 <- numeric(len + 1)
    Pt1 <- numeric(len + 1)
    ## Backward Simulation Container:
    xBS <- numeric(len)
    mut <- numeric(len)
    Lt  <- numeric(len)
    ## Likekihood Estimation Container:
    yt1 <- numeric(len)
    Ft1 <- numeric(len)
    # I. FORWARD FILTERING:
    # 0. Initialization
    x00 <- stateInit
    P00 <- varStateEvol / (1 - phi^2)
    xt1[1] <- phi * x00
    Pt1[1] <- phi^2 * P00 + varStateEvol
    # 1. Iteration: t=1, ...,len
    for (t in 1:len) {
      yt1[t] <- xt1[t]
      Ft1[t] <- Pt1[t] + varObs
      Ft1Inv <- Ft1[t]^(-1)

      xt[t] <- xt1[t] + Pt1[t] * Ft1Inv * (y[t] - yt1[t])
      Pt[t] <- Pt1[t] - Pt1[t] * Ft1Inv * Pt1[t]

      xt1[t + 1] <- phi * xt[t]
      Pt1[t + 1] <- phi^2 * Pt[t] + varStateEvol
    }
    # II. BACKWARD SMOOTHING/SIMULATION
    xBS[len] <- rnorm(1, xt[len], sqrt(Pt[len]))
    for (t in (len - 1):1) {
        tmpVCM <- Pt[t] * phi * (phi^2 * Pt[t] + varStateEvol)^(-1)
        mut[t] <- xt[t] + tmpVCM * (xBS[t + 1] - phi * xt[t])
        Lt[t]  <- Pt[t] - tmpVCM * phi * Pt[t]

        xBS[t] <- rnorm(1, mut[t], sqrt(Lt[t]))
    }
    # III. LOG-LIKELIHOOD COMPUTATION
    llOut1 <- -len / 2 * log(2 * pi)
    llOut2 <- sum(log(Ft1)) + sum((y - yt1)^2 / Ft1)
    llEst  <- llOut1 - 0.5 * llOut2

    return(list(logLikeliEstim = llEst,
                xBackwardSimul = xBS))
}
