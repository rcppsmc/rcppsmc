compareNCestimates <- function(dataY,
                               trueStates = NULL,
                               numParticles = 1000L,
                               simNumber = 100L,
                               modelParameters = list(stateInit = 0,
                                                      phi = 0.7,
                                                      varStateEvol = 1,
                                                      varObs = 1),
                               plot = FALSE) {
    # if no data supplied, use default
    if (missing(dataY)) {
        dataSim  <- simGaussianSSM()
        dataY    <- dataSim$measurements
        # for simulated data, the true states are available
        trueStates <- dataSim$states
    } else if (is.data.frame(dataY) && length(dataY) == 1 ||
          is.matrix(dataY) && ncol(dataY) == 1) {
        dataY <- as.vector(dataY)
    } else if (!is.vector(dataY)) {
        stop("Wrong data format: either dataframe/matrix with single column or vector.")
    }

    resFFBS <- kalmanFFBS(dataY,
                          stateInit = modelParameters$stateInit,
                          phi = modelParameters$phi,
                          varStateEvol = modelParameters$varStateEvol,
                          varObs = modelParameters$varObs)
    resSMC <- compareNCestimates_imp(dataY,
                                     numParticles,
                                     simNum = simNumber,
                                     parInits = modelParameters,
                                     resFFBS$xBackwardSimul)
    if (isTRUE(plot)) {
      par(mfrow = c(3, 2))
      if (!is.null(trueStates)) {
          layoutMatrix <- matrix(c(1, 2,
                                   3, 4), nrow = 2, ncol = 2,
                                   byrow = TRUE)

          layout(mat = layoutMatrix,
                 heights = c(1, 2),   # Heights of the rows
                 widths = c(1, 1))       # Widths of the columns

        title1 <- "measurements (green) vs. true latent states (red)"
        title2 <- "true latent states (red) vs. backward simulation path (blue)"
        plot(dataY, type = "l", col = "forestgreen",
             main = title1, ylab = "", xlab = "time index t")
        lines(trueStates, type = "l", col = "red")
        plot(trueStates, type = "l", col = "red",
             main = title2, ylab = "", xlab = "time index t")
        lines(resFFBS$xBackwardSimul, type = "l", col = "blue")
      } else {
        title3 <- "measurements (green) vs. backward simulation path (blue)"
        plot(dataY, type = "l", col = "forestgreen",
             main = title3, ylab = "", xlab = "time index t")
        lines(resFFBS$xBackwardSimul, type = "l", col = "blue")
        plot.new()
      }
      dataBoxplotsSMC  <- data.frame(llEstimates = as.vector(resSMC$smcOut),
                                     type = c(rep("multinomial", times = simNumber),
                                              rep("residual", times = simNumber),
                                              rep("stratified", times = simNumber),
                                              rep("systematic", times = simNumber)))
      dataBoxplotsCSMC <- data.frame(llEstimates = as.vector(resSMC$csmcOut),
                                     type = c(rep("multinomial", times = simNumber),
                                              rep("residual", times = simNumber),
                                              rep("stratified", times = simNumber),
                                              rep("systematic", times = simNumber)))
      boxplot(llEstimates ~ type, data = dataBoxplotsSMC, col = "bisque",
              ylab = "resampling type",
              xlab = "log-likelihood value",
              main = "Standard BPF log-likelihood estimates",
              sub = "(red: Kalman log-likelihood estimate)",
              horizontal = TRUE)
      abline(v = resFFBS$logLikeliEstim, col = "red")
      boxplot(llEstimates ~ type, data = dataBoxplotsCSMC, col = "bisque",
              ylab = "resampling type",
              xlab = "log-likelihood value",
              main = "Conditional BPF log-likelihood estimates",
              sub = "(red: Kalman log-likelihood estimate)",
              horizontal = TRUE)
      abline(v = resFFBS$logLikeliEstim, col = "red")
    }
    return(list(outSMC = resSMC, outKalman = resFFBS))
}
kalmanFFBS <- function(dataY,
                       stateInit,
                       phi,
                       varStateEvol,
                       varObs) {
    y   <- dataY
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
