## example 2 of Johansen (2009)

rareEventsEx <- function(number=100, iterations=10, threshold=5.0, schedule=30.0) {

    res <- .Call("rareEvents", number, iterations, threshold, schedule, package="RcppSMC")

    invisible(res)
}

