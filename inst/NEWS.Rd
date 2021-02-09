\name{NEWS}
\title{News for Package \pkg{RcppSMC}}
\newcommand{\ghpr}{\href{https://github.com/rcppsmc/rcppsmc/pull/#1}{##1}}
\newcommand{\ghit}{\href{https://github.com/rcppsmc/rcppsmc/issues/#1}{##1}}

\section{Changes in RcppSMC version 0.2.2 (2020-08-30)}{
  \itemize{
    \item Package helper files \code{.editorconfig} added (Adam in \ghpr{43}).
    \item Change const correctness and add return (Leah in \ghpr{44}).
    \item Updates to continuous integration and R versions used (Dirk)
    \item Accomodate CRAN request, other updates to CRAN Policy (Dirk in
    \ghpr{49} fixing \ghit{48}).
  }
}

\section{Changes in RcppSMC version 0.2.1 (2018-03-18)}{
  \itemize{
    \item The sampler now has a copy constructor and assignment overload
    (Brian Ni in \ghpr{28}).
    \item The SMC library component can now be used in header-only mode
    (Martin Lysy in \ghpr{29}).
    \item Plugin support was added for use via \code{cppFunction()} and other
    Rcpp Attributes (or \pkg{inline} functions (Dirk \ghpr{30}).
    \item The sampler copy ctor/assigment operator is now
    copy-constructor safe (Martin Lysy In \ghpr{32}).
    \item A bug in state variance calculation was corrected (Adam in
    \ghpr{36} addressing \ghit{34}).
    \item History getter methods are now more user-friendly (Tiberiu
    Lepadatu in \ghpr{37}).
    \item Use of \code{pow} with atomic types was disambiguated to
    \code{std::pow}) to help the Solaris compiler (Dirk in \ghpr{42}).
  }
}

\section{Changes in RcppSMC version 0.2.0 (2017-08-28)}{
  \itemize{
    \item Also use \code{.registration=TRUE} in \code{useDynLib} in
    \code{NAMESPACE}
    \item Multiple Sequential Monte Carlo extensions (Leah South as part of
    Google Summer of Code 2017)
    \itemize{
      \item Switching to population level objects (\ghpr{2} and \ghpr{3}).
      \item Using Rcpp attributes (\ghpr{2}).
      \item Using automatic RNGscope (\ghpr{4} and \ghpr{5}).
      \item Adding multiple normalising constant estimators (\ghpr{7}).
      \item Static Bayesian model example: linear regression  (\ghpr{10}
      addressing \ghit{9}).
      \item Adding a PMMH example (\ghpr{13} addressing \ghit{11}).
      \item Framework for additional algorithm parameters and adaptation
      (\ghpr{19} addressing \ghit{16}; also \ghpr{24} addressing \ghit{23}).
      \item Common adaptation methods for static Bayesian models
      (\ghpr{20} addressing \ghit{17}).
      \item Supporting MCMC repeated runs (\ghpr{21}).
      \item Adding adaptation to linear regression example (\ghpr{22}
      addressing \ghit{18}).
    }
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
