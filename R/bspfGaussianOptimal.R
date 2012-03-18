
blockpfGaussianOpt <- function(data=c(), particles=1000, lag=5, plot=FALSE) {

    if (length(data == 0)) {
       #Include some error handling here
       return;        
    }
    res <- .Call("blockpfGaussianOpt", data, particles, lag, package="RcppSMC")

    if (plot) {
        time   = 1:length(data);
        mvect  = t(res$weight) %*% res$values / sum(res$weight);
	sqvect = t(res$weight) %*% res$values^2 / sum(res$weight);
	sdvect = sqrt(sqvect - mvect^2);
 
       plot(time, mvect, col='dark red', 'l', lty = 1, lwd=3, xlab = 'Iteration',
       ylab='State', main='Mean and 1, 2 standard deviation credible
       intervals with observations', xlim = c(0,length(data)), ylim=c(min(mvect - 2.1
       * (sdvect)), max(mvect+2.1*sdvect))
       )
       lines(time, mvect + sdvect, lty=3, col='dark blue')
       lines(time, mvect - sdvect, lty=3, col='dark blue')
       lines(time, mvect + 2 * sdvect, lty=2, col='dark blue')
       lines(time, mvect - 2 * sdvect, lty=2, col='dark blue')
       points(time, data, col = 'dark green', cex=0.5)
    }

    invisible(res)
}



