LinRegLA<- function(model, particles=1000,temperatures = seq(0,1,0.05)^5) {

    if (model ==1){
        Data <- cbind(RcppSMC::radiata$y,RcppSMC::radiata$x1)
    } else if (model == 2){
        Data <- cbind(RcppSMC::radiata$y,RcppSMC::radiata$x2)
    } else{
        stop("Please choose a valid model (1 or 2).")
    }

    res <- LinRegLA_impl(as.matrix(Data), as.matrix(temperatures), particles)

    invisible(res)
}
