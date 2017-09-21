LinRegLA_adapt<- function(model, particles = 1000, resampTol = 0.5, tempTol = 0.9) {

    if (model ==1){
        Data <- cbind(RcppSMC::radiata$y,RcppSMC::radiata$x1)
    } else if (model == 2){
        Data <- cbind(RcppSMC::radiata$y,RcppSMC::radiata$x2)
    } else{
        stop("Please choose a valid model (1 or 2).")
    }
    
    res <- LinRegLA_adapt_impl(as.matrix(Data), particles, resampTol, tempTol)

    invisible(res)
}
