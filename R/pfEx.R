
pfEx<- function(filename="", particles=1000, plot=FALSE) {

    if (filename=="") {
        filename <- system.file("sampleData", "pf-data.csv", package="RcppSMC")
    }
    res <- .Call("pf", filename, particles, package="RcppSMC")

    if (plot) {
      ## plot 5.1 from vignette / paper
      data <- read.table(file=filename, skip=1, header=FALSE,
                         col.names=c("x","y"), sep="")
      with(data, plot(x,y,col="red"))
      with(res, lines(Xm, Ym, lty="dashed"))
    }

    invisible(res)
}



