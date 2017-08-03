nonLinPMMH<- function(data, particles=5000, iterations=10000, burnin = 0, plot=FALSE) {

    if (missing(data)) {
         warning("data argument contained no data, using data simulated from the model.")
         data <- simNonlin(len=500,var_init=5,var_evol=10,var_obs=1,cosSeqOffset=0)$data
    }
    
    if (burnin>=iterations) {
        stop("Burn-in must be less than iterations.")
    }
    
    res <- nonLinPMMH_impl(as.matrix(data), particles, iterations)
    
	res.plot <- res[burnin+1:iterations,]
	
    # Replicating Figure 4(c) from Andrieu et al. (2010).
    # May want to add autocorrelation plots of the parameters, like in Figures 5(b) and 5(d).
    if (plot) {
        par(mfrow=c(2,3),oma=c(0,0,2,0))
        with(res.plot, hist(samples_sigv,
                       xlab=expression(sigma_v),ylab='density',
                       main = NA))
        abline(v=sqrt(10),lty=2, col='darkblue',lwd=2)
        with(res.plot, plot(samples_sigw,samples_sigv,
                       xlab=expression(sigma_w),ylab=expression(sigma_v),
                       main = NA, col='darkblue'))
        with(res.plot, plot(samples_sigv,type="l",
                       ylab=expression(sigma_v),
                       main = NA,xlim=c(0,iterations-burnin)))
        with(res.plot, plot(samples_sigv,samples_sigw,
                       xlab=expression(sigma_v),ylab=expression(sigma_w),
                       main = NA, col='darkblue'))
        with(res.plot, hist(samples_sigw,
                       xlab=expression(sigma_w),ylab='density',
                       main = NA))
        abline(v=1,lty=2, col='darkblue',lwd=2)
        with(res.plot, plot(samples_sigw,type="l",
                       ylab=expression(sigma_w),
                       main = NA,xlim=c(0,iterations-burnin)))
        title("Posterior Estimates",outer=TRUE)
        par(mfrow=c(1,1))
    }  

    invisible(res)
}
