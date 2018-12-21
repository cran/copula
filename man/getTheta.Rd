\name{getTheta}
\title{Get the Parameter(s) of a Copula}
\alias{getTheta}
\alias{getTheta-methods}
\alias{getTheta,parCopula-method}
\alias{getTheta,copula-method}
\alias{getTheta,acopula-method}
\alias{getTheta,xcopula-method}
\alias{getTheta,khoudrajiCopula-method}
\alias{getTheta,mixCopula-method}
\alias{getTheta,rotCopula-method}
\description{
  Get the parameter (vector) \eqn{\theta}{theta} (\code{theta}) of
  a copula, see \code{\link{setTheta}} for more background.
}
\usage{%% S4 generic and all S4 methods have same signature (no '...')
getTheta(copula, freeOnly = TRUE, attr = FALSE, named = attr)
}
\arguments{
  \item{copula}{an \R object of class \code{"\linkS4class{parCopula}"},
    e.g., \code{"\linkS4class{copula}"}.}% only, for now
  %%, i.e., any copula from package \pkg{copula}.
  \item{freeOnly}{logical indicating that only non-fixed aka
    \dQuote{free} parameters are to be returned (as vector).}
  \item{attr}{logical indicating if \code{\link{attributes}} (such as
    lower and uppder bounds for each parameters) are to be returned as well.}
  \item{named}{logical if the resulting parameter vector should have
    \code{\link{names}}.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(copula = "parCopula")}}{a default method,
      returning \code{numeric(0).}}
    \item{\code{signature(copula = "copula")}}{ .. }
    \item{\code{signature(copula = "acopula")}}{ .. }
    \item{\code{signature(copula = "khoudrajiCopula")}}{ .. }
    \item{\code{signature(copula = "mixCopula")}}{ .. }
    \item{\code{signature(copula = "rotCopula")}}{ .. }
    \item{\code{signature(copula = "xcopula")}}{ .. }
  }
}
%% \details{
%% }
\value{
  parameter vector of the copula, a \code{\link{numeric}} vector,
  possibly with names and other attributes (depending on the \code{attr}
  and \code{named} arguments).
}
\seealso{
  \code{\link{setTheta}}, its inverse.
}
\examples{
showMethods("getTheta")
\dontshow{stopifnot(all.equal(}
getTheta(setTheta(copClayton, 0.5)) # is 0.5
\dontshow{, 0.5))}
}
\keyword{manip}
\keyword{methods}
