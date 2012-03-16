pfNonlinBS <- function(data=c(), particles=500, plot=FALSE) {
    if (length(data == 0)) {
       #Include some error handling here
       return;        
    }
    res <- .Call("pfNonlinBS", data, particles, package="RcppSMC")

    time <- 1:length(data);
    if (plot) {
      with(res, plot(time,mean,type='l',xlab='time',ylab='estimate'))
      with(res,polygon(c(time,51-time),c(mean-2*sd,(mean+2*sd)[seq(length(data),1,-1)]),col='palegreen1',border=NA))
      with(res,polygon(c(time,51-time),c(mean-1*sd,(mean+1*sd)[seq(length(data),1,-1)]),col='palegreen3',border=NA))
      with(res,lines(time,mean, lwd=2, col='darkblue'))
    }

    invisible(res)
}



