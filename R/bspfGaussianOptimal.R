#' Block Sampling Particle Filter (Linear Gaussian Model; Optimal Proposal).
#'
#' @description The \code{blockpfGaussianOpt} function provides a simple example for
#' \pkg{RcppSMC}. It is based on a block sampling particle filter for a linear
#' Gaussian model. This is intended only to illustrate the potential of block
#' sampling; one would not ordinarily use a particle filter for a model in
#' which analytic solutions are available. The 'optimal' block sampler in the
#' sense of Doucet, Briers and Senecal (2006) can be implemented in this case.
#' 
#' @param data	A vector variable containing the sequence of observations.
#' 
#' @param particles	An integer specifying the number of particles.
#' 
#' @param lag An integer specifying the length of block to use.
#' 
#' @param plot A boolean variable describing whether plot should
#' illustrate the estimated path along with the uncertainty.
#' 
#' @return	The \code{blockpfGaussianOpt} function returns a matrix containing the final
#' sample paths and a vector containing their weights. The logarithm of the
#' estimated ratio of normalising constants between the final and initial
#' distributions is also returned.
#'
#' @details The \code{blockpfGaussianOpt} function provides a simple example for
#' \pkg{RcppSMC}. It is based on a simple linear Gaussian state space model in
#' which the state evolution and observation equations are:
#' 	x(n) = x(n-1) + e(n) and 
#' 	y(n) = x(n) + f(n)
#' where e(n) and f(n) are mutually-independent standard normal random
#' variables. The 'optimal' block-sampling proposal described by Doucet
#' et al (2006) is employed. 
#' 
#' @references	A. Doucet, M. Briers, and S. Senecal. Efficient Block Sampling Strategies
#' for sequential Monte Carlo methods. Journal of Computational and Graphical
#' Statistics, 15(3):693-711, 2006.
#' 
#' @examples
#' sim <- simGaussian(len=250)
#' res <- blockpfGaussianOpt(sim$data,lag=5,plot=TRUE)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @export
blockpfGaussianOpt <- function(data=c(), particles=1000, lag=5, plot=FALSE) {

    if (length(data) == 0) {
       warning("data argument contained no data, using data simulated from the model.")
       data <- simGaussian(len=250)$data
    }
    res <- blockpfGaussianOpt_impl(data, particles, lag)

    if (plot) {
        time   = 1:length(data);
        mvect  = t(res$weight) %*% res$values / sum(res$weight);
	sqvect = t(res$weight) %*% res$values^2 / sum(res$weight);
	sdvect = sqrt(sqvect - mvect^2);

       plot(time, mvect, 'l', lty = 1, lwd=3, xlab = 'Iteration', ylab='State',
            main='Mean and 1, 2 standard deviation credible intervals with observations',
	    xlim = c(0,length(data)), ylim=c(min(mvect - 2.1 * (sdvect)), max(mvect+2.1*sdvect)))

       polygon(c(time,seq(length(data),1,-1)),c(mvect-2*sdvect,(mvect+2*sdvect)[seq(length(data),1,-1)]),col='palegreen1',border=NA)
       polygon(c(time,seq(length(data),1,-1)),c(mvect-1*sdvect,(mvect+1*sdvect)[seq(length(data),1,-1)]),col='palegreen3',border=NA)
       lines(time, mvect, lwd=2, col='dark blue')

       points(time, data, col = 'dark green', cex=0.5)
    }

    invisible(res)
}



