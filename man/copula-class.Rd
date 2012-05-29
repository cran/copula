\name{copula-class}
\title{Class "copula"}
\docType{class}
\alias{copula-class}
\alias{Copula-class}
\alias{plackettCopula-class}
\alias{dim,copula-method}
\description{A class representing multivariate distributions with uniform margins.}
\section{Objects from the Class}{
  Objects could be created by calls of the form \code{new("copula", ...)},
  but typically are by \code{\link{tCopula}()},
  \code{\link{evCopula}()}, etc.

  Note that there's also a virtual class \code{"Copula"}, the union (see
  \code{\link{setClassUnion}}) of \code{"copula"} and class
  \code{"\linkS4class{nacopula}"}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{the dimension of the copula, of class \code{"integer"}.}
    \item{\code{parameters}:}{numeric parameter vector (class \code{"numeric"}).}
    \item{\code{param.names}:}{\code{"character"} vector of parameter names.}
    \item{\code{param.lowbnd}:}{lower bounds for the parameters, of
      class \code{"numeric"}.}
    \item{\code{param.upbnd}:}{upper bounds for the parameters, of class
      \code{"numeric"}.}
    \item{\code{message}:}{Object of class \code{"character"}, family names
      of the copula.}
  }
}
\note{
  The \code{"copula"} class is extended by the
  \code{"\linkS4class{evCopula}"}, \code{"\linkS4class{archmCopula}"},
  and \code{"\linkS4class{ellipCopula}"} classes.  Instances of
  such copulas can be created via functions \code{\link{evCopula}},
  \code{\link{archmCopula}} and \code{\link{ellipCopula}}.

  \code{"plackettCopula"} and \code{"\linkS4class{fgmCopula}"} are special types of
  copulas which do not belong to either one of the three classes above.
}
\section{Warning}{
  This implementation is still at the experimental
  stage and is subject to change during the development.
}
\seealso{
  Help for the (sub)classes
  \code{\linkS4class{archmCopula}},
  \code{\linkS4class{ellipCopula}},
  \code{\linkS4class{evCopula}}, and
  \code{\linkS4class{fgmCopula}}.

  The Archimedean and nested Archimedean classes (from former package
  \pkg{nacopula}), with a more extensive list of slots (partly instead
  of methods),
  \code{\linkS4class{acopula}}, and \code{\linkS4class{nacopula}}.
}
\examples{
hc <- evCopula("husler", 1.25)
dim(hc)
smoothScatter(u <- rcopula(hc, 2^11))
tailIndex   (hc)
kendallsTau (hc)
spearmansRho(hc)
str(hc)
}
\keyword{classes}
