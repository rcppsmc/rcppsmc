
pfEx<- function(data, particles=1000, plot=FALSE, onlinePlot) {

    # if no data supplied, use default
    if (missing(data)) data <- getEx1Data()

    if (missing(onlinePlot)) {
        useOnline <- FALSE
        onlinePlot <- function() { NULL }
    } else {
        useOnline <- TRUE
        # set up x11
        x11(width=3,height=3)
        par(mar=c(3,3,1,1),cex=0.8, pch=19)
    }

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"),
              class(onlinePlot) == "function")

    res <- .Call("pf", as.matrix(data),
                 particles,
                 useOnline,
                 onlinePlot,
                 package="RcppSMC")

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

pfExOnlinePlot <- function(xm, ym) {
    plot(xm, ym, ylim=c(-7,0), xlim=c(2,14))
    Sys.sleep(0.05)
}

