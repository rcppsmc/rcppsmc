\name{RcppSMC.package.skeleton}
\alias{RcppSMC.package.skeleton}
\title{
Create a skeleton for a new package that intends to use RcpSMCp
}
\description{
	\code{RcppSMC.package.skeleton} automates the creation of
	a new source package that intends to use features of RcppSMC.

	It is based on the \link[utils]{package.skeleton} and
	\link[pkgKitten]{kitten} (from \link[pkgKitten]{pkgKitten}) functions, the
	latter being a Wrapper around \link[utils]{package.skeleton} to make a
	package pass \samp{R CMD check}	without	complaints. If
	\link[pkgKitten]{pkgKitten} is not installed,
	\link[utils]{package.skeleton} is executed instead.
}
\usage{
RcppSMC.package.skeleton(name = "anRpackage", list = character(),
	environment = .GlobalEnv, path = ".")
}
\arguments{
	\item{name}{See \link[utils]{package.skeleton}}
	\item{list}{See \link[utils]{package.skeleton}}
	\item{environment}{See \link[utils]{package.skeleton}}
	\item{path}{See \link[utils]{package.skeleton}}
}
\details{
	In addition to \link[utils]{package.skeleton} :

	The \samp{DESCRIPTION} file gains a Depends line requesting that
	the package depends on Rcpp, RcppArmadillo and RcppSMC and
	a LinkingTo line so that the package finds the associated header files.

	The \samp{NAMESPACE}, if any, gains a \code{useDynLib} directive.

	The \samp{src} directory is created if it does not exists and
	a \samp{Makevars} file is added setting the environment variable
	\samp{PKG_LIBS} to accomodate the necessary flags
	to link with the Rcpp library.

	An example file \samp{rcppsmc_hello_world.cpp} is created in the \samp{src}. An R file \samp{rcppsmc_hello_world.R} is
	expanded in the \samp{R} directory, the \code{rcppsmc_hello_world} function
	defined in this files makes use of the C++ function \samp{rcppsmc_hello_world}
	defined in the C++ file. These files are given as an example and should
	eventually by removed from the generated package.
}
\value{
Nothing, used for its side effects
}
\seealso{
\link[utils]{package.skeleton}
\link[pkgKitten]{kitten}
}
\references{
Read the \emph{Writing R Extensions} manual for more details.

Once you have created a \emph{source} package you need to install it:
see the \emph{R Installation and Administration} manual,
\code{\link{INSTALL}} and \code{\link{install.packages}}.
}
\examples{
\dontrun{
RcppSMC.package.skeleton( "foobar" )
}
}
\keyword{ programming }
