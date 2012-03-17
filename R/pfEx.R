
pfEx<- function(data, particles=1000, plot=FALSE) {

    # if no data supplied, use default
    if (missing(data)) data <- getEx1Data()

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"))

    res <- .Call("pf", as.matrix(data), particles, package="RcppSMC")

    if (plot) {
      ## plot 5.1 from vignette / paper
      with(data, plot(x, y, col="red"))
      with(res, lines(Xm, Ym, lty="dashed"))
    }

    invisible(res)
}

# simple convenience function, should probably make the data a
# data component of the package...
getEx1Data <- function() {
    file <- system.file("sampleData", "pf-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x","y"))
    invisible(dat)
}


