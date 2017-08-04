LinReg<- function(model, particles=1000, plot=FALSE) {

    data(radiata, envir = environment())
    if (model ==1){
        Data <- cbind(radiata$y,radiata$x1)
    } else if (model == 2){
        Data <- cbind(radiata$y,radiata$x2)
    } else{
        stop("Please choose a valid model (1 or 2).")
    }

    res <- LinReg_impl(as.matrix(Data),particles)

    if (plot) {
        par(mfrow=c(1,3),oma=c(0,0,2,0))
        with(res, plot(density(theta[,1],weights=weights),type='l',
                       xlab=expression(alpha),ylab='density',
                       main = NA, col='darkblue'))
        with(res, plot(density(theta[,2],weights=weights),type='l',
                       xlab=expression(beta),ylab='density',
                       main = NA, col='darkblue'))
        with(res, plot(density(theta[,3],weights=weights),type='l',
                       xlab=expression(phi),ylab='density',
                       main = NA, col='darkblue'))
        title("Posterior Estimates",outer=TRUE)
        par(mfrow=c(1,1))
    }

    invisible(res)
}

## silence a NOTE from 'R CMD check --as-cran'
utils::globalVariables(c("radiata"))
