#' Particle Filter Example.
#'
#' @description The \code{pfLineartBS} function provides a simple example for
#' \pkg{RcppSMC}. It is based on the first example in \code{SMCTC} and
#' the discussion in Section 5.1 of Johansen (2009). A simple 'vehicle
#' tracking' problem of 100 observations is solved with 1000 particles.
#' 
#' @param dat	A two-column matrix or dataframe containing x and y
#' values. The default data set from Johansen (2009) is used as the
#' default if no data is supplied.
#' 
#' @param particles	An integer specifying the number of particles.
#' 
#' @param plot	A boolean variable describing whether plot should
#' illustrate the estimated path along with the data.
#' 
#' @param onlinePlot	A user-supplied callback function which is called with the
#' x and y position vectors during each iteration of the algorithm; see
#' pfExOnlinePlot for a simple example.
#' 
#' @return	The \code{pfLineartBS} function returns a \code{data.frame} containing as many rows as in
#' the input data, and four columns corresponding to the estimated \eqn{x}{x} and
#' \eqn{y}{y} coordinates as well as the estimated velocity in these two
#' directions.
#'
#' @details The \code{pfLineartBS} function provides a simple example for
#' \pkg{RcppSMC}. The model is linear with t-distributed innovations.
#' It is based on the \code{pf} example in the
#' \code{SMCTC} library, and discussed in the Section 5.1 of his
#' corresponding paper (Johansen, 2009).
#' 
#' @references	A. M. Johansen. SMCTC: Sequential Monte Carlo in C++.
#' Journal of Statistical Software, 30(6):1-41, April
#' 2009. \url{http://www.jstatsoft.org/v30/i06/paper}
#' 
#' @seealso The SMCTC paper and code at \url{http://www.jstatsoft.org/v30/i06/paper}.
#' 
#' @examples
#' res <- pfLineartBS(plot=TRUE)
#' if (interactive()) ## if not running R CMD check etc
#' res <- pfLineartBS(onlinePlot=pfLineartBSOnlinePlot)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @export
pfLineartBS<- function(dat, particles=1000, plot=FALSE, onlinePlot) {

    # if no data supplied, use default
    if (missing(dat)){
		data(pfdata)
		dat <- pfdata
	}

    if (missing(onlinePlot)) {
        useOnline <- FALSE
        onlinePlot <- function() { NULL }
    } else {
        useOnline <- TRUE
        # set up graphics window
        dev.new(width=3,height=3)
        par(mar=c(3,3,1,1),cex=0.8, pch=19, ask=FALSE)
    }

    # more eloquent tests can be added
    stopifnot(nrow(dat) > 0,
              ncol(dat) == 2,
              colnames(dat) == c("x", "y"),
              class(onlinePlot) == "function")

    res <- pfLineartBS_impl(as.matrix(dat),
                 particles,
                 useOnline,
                 onlinePlot)

    if (plot) {
        ## plot 5.1 from vignette / paper
        with(dat, plot(x, y, col="red"))
        with(res, lines(Xm, Ym, lty="dashed"))
    }

    invisible(res)
}
pfLineartRange <- function(rrng)
{
   min <- rrng[1]
   max <- rrng[2]

   if(min > 0) {
      rmin = exp(floor(log(min)))
   } else if (min < 0) {
      rmin = -exp(ceiling(log(-min)))
   } else {
      rmin = 0
   }

   if(max > 0) {
      rmax = exp(ceiling(log(max)))
   } else if (max < 0){
      rmax = exp(floor(log(-max)))
   } else {
      rmax = 0;
   }

   invisible(c(rmin,rmax))
}

#' For online plotting example
#'
#' The \code{pfLineartBSOnlinePlot} function provides a simple default
#' \sQuote{online} plotting function that is invoked during the
#' estimation process.
#' 
#' @param xm Vector with x position.
#' @param ym Vector with y position.
#' 
#' @details Using the simple \code{pfExOnlinePlot} function illustrates how
#' callbacks into R, for example for plotting, can be made during the
#' operation of SMC algorithm.
#' 
#' @rdname pfLineartBS
#' @export
pfLineartBSOnlinePlot <- function(xm, ym) {
    plot(xm, ym, xlim = pfLineartRange(range(xm)), ylim=pfLineartRange(range(ym)))
    # FIXME sleep time should also be a variable
    Sys.sleep(0.05)
}

