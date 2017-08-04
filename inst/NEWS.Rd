\name{NEWS}
\title{News for Package \pkg{RcppSMC}}
\newcommand{\ghpr}{\href{https://github.com/eddelbuettel/rcppsmc/pull/#1}{##1}} 
\newcommand{\ghit}{\href{https://github.com/eddelbuettel/rcppsmc/issues/#1}{##1}} 

\section{Changes in RcppSMC version 0.1.6 [unreleased] (2017-xx-yy)}{
  \itemize{
    \item Also use \code{.registration=TRUE} in \code{useDynLib} in
    \code{NAMESPACE}
    \item The particles are now represented internally as a population 
    level object rather than particle level objects. As a result, 
    the user defined initialisation, move and MCMC functions 
    take a different format (Leah South in \ghpr{2} and \ghpr{3}). 
    \item Rcpp attributes are now used (Leah South in \ghpr{2}). 
    \item The automatic RNGscope from Rcpp attributes is used 
    (Leah South in \ghpr{4} and \ghpr{5}). 
    \item Estimation of the log ratio of normalising constants 
    can now be performed via the standard SMC estimator and 
    multiple path sampling methods (Leah South in \ghpr{7}) 
    \item A linear regression example has been added to  
    demonstrate how the package can be used to perform  
    inference for static Bayesian models and 
    Bayesian model choice (Leah South in \ghit{9} 
    and \ghpr{10}). 
    \item An example of particle marginal Metropolis
    Hastings has been added (Leah South in \ghit{11} 
    and \ghpr{13}). 
  }
}
  
\section{Changes in RcppSMC version 0.1.5 (2017-03-03)}{
  \itemize{
    \item Correct \code{.Call} to use \code{PACKAGE=} argument
    \item \code{DESCRIPTION}, \code{NAMESPACE}, \code{README.md} changes
    to comply with current \code{R CMD check} levels
    \item Added file \code{init.c} with calls to \code{R_registerRoutines()}
    \code{and R_useDynamicSymbols()}
    \item Updated \code{.travis.yml} file for continuous integration
  }
}

\section{Changes in RcppSMC version 0.1.4 (2014-05-02)}{
  \itemize{
    \item Added missing delete operator to destructor in \code{sampler} 
  }
}

\section{Changes in RcppSMC version 0.1.3 (2014-05-01)}{
  \itemize{
    \item Bugfix in \code{Sampler.iterate()} for memory overrun detected by valgrind
  }
}

\section{Changes in RcppSMC version 0.1.2 (2014-02-06)}{
  \itemize{
    \item Updated for \CRANpkg{Rcpp} 0.11.0 with explicit \code{importFrom}
    in \code{NAMESPACE} and corresponding versioned \code{Imports:} in \code{DESCRIPTION};
    also removed linking instruction from \code{src/Makevars} as no it
    is longer needed with this new \CRANpkg{Rcpp} version
    \item Addded GitHub / Travis CI support
    \item Use more portable \code{dev.new()} rather than \code{x11()} in \code{pfLinearBS.R}
    \item Applied some corrections to \code{pfNonlinBS.R} example
    \item Converted \code{NEWS} to \code{NEWS.Rd}
  }
}

\section{Changes in RcppSMC version 0.1.1 (2012-05-14)}{
  \itemize{
    \item Minor \code{g++-4.7} build fix of using \code{std::max()} explicitly
  }
}

\section{Changes in RcppSMC version 0.1.0 (2012-03-21)}{
  \itemize{
    \item Version 0.1.0 and initial CRAN upload
  }
}
