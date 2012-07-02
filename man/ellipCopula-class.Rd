\name{ellipCopula-class}
\docType{class}
\alias{ellipCopula-class}
\alias{normalCopula-class}
\alias{tCopula-class}
\title{Class "ellipCopula"}
\description{Copulas generated from elliptical multivariate distributions. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ellipCopula", ...)},
  or by function \code{\link{ellipCopula}}.
}
\section{Slots}{
  \describe{
    \item{\code{dispstr}:}{Object of class \code{"character"}, indicating
    how the dispersion matrix is parameterized; can 'ex', 'ar1', 'toep',
    or 'un'.}
    \item{\code{dimension}:}{Object of class \code{"numeric"}, dimension
      of the copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"},
      parameter value.}
    \item{\code{param.names}:}{Object of class \code{"character"},
      parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"},
      parameter lower bound. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"},
      parameter upper bound. }
    \item{\code{fullname}:}{Object of class \code{"character"}, family names
      of the copula.}
  }
}
\section{Extends}{
  Class \code{"ellipCopula"} extends class \code{"\linkS4class{copula}"}
  directly.Class \code{"normalCopula"} and \code{"tCopula"} extends class
  \code{"ellipCopula"} directly.
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "tCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "tCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "tCopula")}: ... }
  }
}
\seealso{
  \code{\link{ellipCopula}},
  \code{\link{copula-class}}.
}
\keyword{classes}
