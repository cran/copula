\name{ellipCopula-class}
\title{Class "ellipCopula" of Elliptical Copulas}
\docType{class}
\alias{ellipCopula-class}
\alias{normalCopula-class}
\alias{tCopula-class}
\description{
  Copulas generated from elliptical multivariate distributions, notably
  Normal- and t-copulas (of specific class \code{"normalCopula"} or
  \code{"tCopula"}, respectively).
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{ellipCopula}()},
   \code{\link{normalCopula}()}, or \code{\link{tCopula}()}.
}
\section{Slots}{
  \describe{
    \item{\code{dispstr}:}{\code{"\link{character}"} string indicating
      how the dispersion matrix is parameterized; one of \code{"ex"},
      \code{"ar1"}, \code{"toep"}, or \code{"un"}, see the \code{dispstr}
      argument of \code{\link{ellipCopula}()}.}
    \item{\code{dimension}:}{Object of class \code{"numeric"}, dimension
      of the copula. }
    \item{\code{parameters}:}{a \code{\link{numeric}}, (vector of) the parameter
      value(s).}
    \item{\code{param.names}:}{\code{\link{character}} vector with names
      for the \code{parameters} slot, of the same length.}
    \item{\code{param.lowbnd}:}{\code{\link{numeric}} vector of lower
      bounds for the \code{parameters} slot, of the same length.}
    \item{\code{param.upbnd}:}{upper bounds for \code{parameters},
      analogous to \code{parm.lowbnd}.}
    \item{\code{fullname}:}{\bold{deprecated}; object of class \code{"character"}, family names
      of the copula.}
  }
}
\section{Extends}{
  Class \code{"ellipCopula"} extends class \code{"\linkS4class{copula}"}
  directly. Classes \code{"normalCopula"} and \code{"tCopula"} extend
  \code{"ellipCopula"} directly.
}
\section{Methods}{
  Many methods are available, notably \code{\link{dCopula}},
  \code{\link{pCopula}}, and \code{\link{rCopula}}.
  Use, e.g., \code{\link{methods}(class = "tCopula")} to find others.
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{
  \code{\link{ellipCopula}} which also documents \code{tCopula()} and
  \code{\link{normalCopula}()};
  \code{\link{copula-class}}.
}
\keyword{classes}
