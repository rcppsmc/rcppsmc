
pfLineartBS<- function(data, particles=1000, plot=FALSE, onlinePlot) {

    # if no data supplied, use default
    if (missing(data)) data <- getPfLineartBSData()

    if (missing(onlinePlot)) {
        useOnline <- FALSE
        onlinePlot <- function() { NULL }
    } else {
        useOnline <- TRUE
        # set up x11
        x11(width=3,height=3)
        par(mar=c(3,3,1,1),cex=0.8, pch=19, ask=FALSE)
    }

    # more eloquent tests can be added
    stopifnot(nrow(data) > 0,
              ncol(data) == 2,
              colnames(data) == c("x", "y"),
              class(onlinePlot) == "function")

    res <- .Call("pfLineartBS", as.matrix(data),
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
getPfLineartBSData <- function() {
    file <- system.file("sampleData", "pf-data.csv", package="RcppSMC")
    dat <- read.table(file, skip=1, header=FALSE, col.names=c("x","y"))
    invisible(dat)
}

pfLineartBSOnlinePlot <- function(xm, ym) {
    # FIXME: xlim and ylim should be dependent on data, but of course
    # the online algorithm does not "know" the full dataset as it
    # works its way through
    plot(xm, ym, xlim=c(-7,0), ylim=c(2,14))
    # FIXME sleep time should also be a variable
    Sys.sleep(0.05)
}

