#' Nonlinear Bootstrap Particle Filter (Univariate Non-Linear State Space Model).
#'
#' @description The \code{pfNonlinBS} function provides a simple example for
#' \pkg{RcppSMC}. It is a simple \dQuote{bootstrap} particle filter which employs
#' multinomial resampling after each iteration applied to the ubiquitous "nonlinear
#' state space model" following Gordon, Salmond and Smith (1993).
#' 
#' @param data	A vector variable containing the sequence of observations.
#' 
#' @param particles	An integer specifying the number of particles.
#' 
#' @param plot	A boolean variable describing whether a plot should
#' illustrate the (posterior mean) estimated path along with one and two
#' standard deviation intervals.
#' 
#' @return  The \code{pfNonlinBS} function returns two vectors, the first containing the posterior
#' filtering means; the second the posterior filtering standard deviations.
#'
#' @details The \code{pfNonlinbs} function provides a simple example for
#' \pkg{RcppSMC}. It is based on a simple nonlinear state space model in
#' which the state evolution and observation equations are:
#' 	x(n) = 0.5 x(n-1) + 25 x(n-1) / (1+x(n-1)^2) + 8 cos(1.2(n-1))+ e(n) and 
#'  y(n) = x(n)^2 / 20 + f(n)
#' where e(n) and f(n) are mutually-independent normal random
#' variables of variances 10.0 and 1.0, respectively. A boostrap proposal
#' (i.e. sampling from the state equation) is used, together with multinomial
#' resampling after each iteration. 
#' 
#' @references	N. J. Gordon, S. J. Salmond, and A. F. M. Smith. Novel approach to
#' nonlinear/non-Gaussian Bayesian state estimation. IEE Proceedings-F, 
#' 140(2):107-113, April 1993.
#' 
#' @examples
#' sim <- simNonlin(len=50)
#' res <- pfNonlinBS(sim$data,particles=500,plot=TRUE)
#' 
#' @author		Adam M. Johansen, Dirk Eddelbuettel and Leah F. South
#' @keywords	programming
#' @export
pfNonlinBS <- function(data, particles=500, plot=FALSE) {
    if (missing(data)) {
         warning("data argument contained no data, using data simulated from the model.")
         data <- simNonlin(len=50)$data
    }
    res <- pfNonlinBS_impl(data,particles)

    time <- 1:length(data);
    if (plot) {
        with(res, plot(time,mean,type='l',
                       main='Filtering Mean and +/- 1,2 standard deviation intervals',
                       xlab='time',ylab='estimate',
                       xlim = c(0,length(data)),
                       ylim = c(min(mean-2.1*sd),max(mean+2.1*sd))))
        with(res,polygon(c(time,seq(length(data),1,-1)),
                         c(mean-2*sd,(mean+2*sd)[seq(length(data),1,-1)]),
                         col='palegreen1',border=NA))
        with(res,polygon(c(time,seq(length(data),1,-1)),
                         c(mean-1*sd,(mean+1*sd)[seq(length(data),1,-1)]),
                         col='palegreen3',border=NA))
        with(res,lines(time,mean, lwd=2, col='darkblue'))
    }

    invisible(res)
}



