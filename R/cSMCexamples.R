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
                          varObs = modelParameters$varObs,
                          simNumber = simNumber)
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
        title2 <- "true latent states (red) vs. avg. backward simulation path (blue)"
        plot(dataY, type = "l", col = "forestgreen",
             main = title1, ylab = "", xlab = "time index t")
        lines(trueStates, type = "l", col = "red")
        plot(trueStates, type = "l", col = "red",
             main = title2, ylab = "", xlab = "time index t")
        lines(rowMeans(resFFBS$xBackwardSimul), type = "l", col = "blue")
      } else {
        title3 <- "measurements (green) vs. avg. backward simulation path (blue)"
        plot(dataY, type = "l", col = "forestgreen",
             main = title3, ylab = "", xlab = "time index t")
        lines(rowMeans(resFFBS$xBackwardSimul), type = "l", col = "blue")
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
                       varObs,
                       simNumber) {
    len <- length(dataY)
    xBS <- matrix(0, nrow = len, ncol = simNumber)
    outKF <- FKF::fkf(a0 = stateInit,
                      P0 = matrix(varStateEvol / (1 - phi^2)),
                      dt = matrix(0, nrow = 1, ncol = 1),
                      ct = matrix(0, nrow = 1, ncol = 1),
                      Tt = array(phi, c(1, 1, 1)),
                      Zt = array(1, c(1, 1, 1)),
                      HHt = array(varStateEvol, c(1, 1, 1)),
                      GGt = array(varObs, c(1, 1, 1)),
                      yt = matrix(dataY, nrow = 1))
    outKS <- FKF::fks(outKF)
    # II. BACKWARD SMOOTHING/SIMULATION
    for (n in 1:simNumber) {
        xBS[len, n] <- rnorm(1, outKS$ahat[len], sqrt(outKS$Vt[, , len]))
        for (t in (len - 1):1) {
            xBS[t, n] <- rnorm(1, outKS$ahat[t], sqrt(outKS$Vt[, , t]))
        }
    }
    # III. LOG-LIKELIHOOD COMPUTATION
    llEst  <- outKF$logLik

    return(list(logLikeliEstim = llEst,
                xBackwardSimul = xBS))
}
