RcppSMC.package.skeleton <- function (name = "anRpackage", list = character(),
                                      environment = .GlobalEnv,path = ".") {
    env <- parent.frame(1)
    if (!length(list)) {
        fake <- TRUE
        assign("Rcpp.fake.fun", function() {
        }, envir = env)
    }
    else {
        fake <- FALSE
    }
    haveKitten <- requireNamespace("pkgKitten", quietly = TRUE)
    skelFunUsed <- ifelse(haveKitten, pkgKitten::kitten, package.skeleton)
    skelFunName <- ifelse(haveKitten, "kitten", "package.skeleton")
    message("\nCalling ", skelFunName, " to create basic package.")
    call <- match.call()
    call[[1]] <- skelFunUsed
    if ("example_code" %in% names(call)) {
        call[["example_code"]] <- NULL
    }
    if (!haveKitten) {
        if (fake) {
            call[["list"]] <- "Rcpp.fake.fun"
        }
    }
    tryCatch(eval(call, envir = env), error  = function(e) {
        cat(paste(e, "\n"))
        stop(paste("error while calling `", skelFunName, "`",
            sep = ""))
    })
    message("\nAdding RcppSMC settings")
    root <- file.path(path, name)
    DESCRIPTION <- file.path(root, "DESCRIPTION")
    if (file.exists(DESCRIPTION)) {
        x <- cbind(read.dcf(DESCRIPTION),
                   Imports = sprintf("Rcpp (>= %s)",
                   packageDescription("Rcpp")[["Version"]]),
                   LinkingTo = "Rcpp, RcppArmadillo, RcppSMC")
        write.dcf(x, file = DESCRIPTION)
        message(" >> added Imports: Rcpp")
        message(" >> added LinkingTo: Rcpp, RcppArmadillo and RcppSMC")
    }
    NAMESPACE <- file.path(root, "NAMESPACE")
    lines <- readLines(NAMESPACE)
    lines <- lines[!grepl("^export.*fake\\.fun", lines)]
    if (!any(grepl("^exportPattern", lines))) {
        lines <- c(lines, "exportPattern(\"^[[:alpha:]]+\")")
    }
    if (!any(grepl("useDynLib", lines))) {
        lines <- c(sprintf("useDynLib(%s, .registration=TRUE)",
            name), "importFrom(Rcpp, evalCpp)", lines)
        writeLines(lines, con = NAMESPACE)
        message(" >> added useDynLib and importFrom directives to NAMESPACE")
    }
    src <- file.path(root, "src")
    if (!file.exists(src)) {
        dir.create(src)
    }
    man <- file.path(root, "man")
    if (!file.exists(man)) {
        dir.create(man)
    }
    skeletonArma <- system.file("skeleton", package = "RcppArmadillo")
    skeletonSMC  <- system.file("skeleton", package = "RcppSMC")
    Makevars <- file.path(src, "Makevars")
    if (!file.exists(Makevars)) {
        file.copy(file.path(skeletonArma, "Makevars"), Makevars)
        message(" >> added Makevars file for baseline compiler flags")
    }
    Makevars.win <- file.path(src, "Makevars.win")
    if (!file.exists(Makevars.win)) {
        file.copy(file.path(skeletonArma, "Makevars.win"), Makevars.win)
        message(" >> added Makevars file for baseline compiler flags (Windows)")
    }

    file.copy(file.path(skeletonSMC, "rcppsmc_hello_world.cpp"), src)
    message(" >> added example src file using RcppSMC functions/classes")

    file.copy(file.path(skeletonSMC, "rcppsmc_hello_world.Rd"), man)
    message(" >> added example Rd file for using RcppSMC functions/classes")

    Rcpp::compileAttributes(root)
    message(" >> invoked Rcpp::compileAttributes to create wrappers")

    if (fake) {
        rm("Rcpp.fake.fun", envir = env)
        unlink(file.path(root, "R", "Rcpp.fake.fun.R"))
        unlink(file.path(root, "man", "Rcpp.fake.fun.Rd"))
    }
    invisible(NULL)
}
