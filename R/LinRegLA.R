LinRegLA<- function(model, particles=1000,temperatures = seq(0,1,0.05)^5) {

    if (model ==1){
        Data <- cbind(radiata$y,radiata$x1)
    } else if (model == 2){
        Data <- cbind(radiata$y,radiata$x2)
    } else{
        stop("Please choose a valid model (1 or 2).")
    }

    res <- LinRegLA_impl(as.matrix(Data), as.matrix(temperatures), particles)

    invisible(res)
}

## silence a NOTE from 'R CMD check --as-cran'
utils::globalVariables(c("radiata"))
