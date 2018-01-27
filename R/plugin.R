## Copyright (C) 2018  Dirk Eddelbuettel
##
## This file is part of RcppSMC.
##
## RcppSMC is free software: you can redistribute it and/or modify it
## under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
##
## RcppSMC is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with RcppSMC.  If not, see <http://www.gnu.org/licenses/>.

inlineCxxPlugin <- function(...) {
    ismacos <- Sys.info()[["sysname"]] == "Darwin"
    openmpflag <- if (ismacos) "" else "$(SHLIB_OPENMP_CFLAGS)"
    plugin <- Rcpp::Rcpp.plugin.maker(include.before = "#include <RcppSMC.h>",
                                      libs           = paste(openmpflag,
                                                             "$(LAPACK_LIBS) $(BLAS_LIBS) $(FLIBS)"),
                                      package        = "RcppSMC")
    settings <- plugin()
    settings$env$PKG_CPPFLAGS <- paste0("-I../inst/include ",
                                       "-I", system.file("include", package="RcppArmadillo"), " ",
                                       openmpflag)
    if (!ismacos) settings$env$USE_CXX11 <- "yes"
    settings
}
