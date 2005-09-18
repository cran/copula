\name{ellipCopula-class}
\docType{class}
\alias{ellipCopula-class}
\alias{normalCopula-class}
\alias{tCopula-class}
\title{Class "ellipCopula" }
\description{ Copulas generated from elliptical multivariate
  distributions. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("ellipCopula",
    ...)}, or by function 'ellipCopula'.
}
\section{Slots}{
  \describe{
    \item{\code{dispstr}:}{Object of class \code{"character"}, indicating
    the type of the dispersion matrix such as 'ex', 'ar1', 'toep', or 'un'.}
    \item{\code{dimension}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.names}:}{Object of class \code{"character"} ~~ }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{message}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Extends}{
  Class \code{"ellipCopula"} extends class \code{"copula"} directly.
  Class \code{"normalCopula"} and \code{"tCopula"} extends class
  \code{"ellipCopula"} directly.
}
\section{Methods}{
  \describe{
    \item{dcopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "normalCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "tCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "tCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "tCopula")}: ... }
  }
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\seealso{
  \code{\link{ellipCopula}}, 
  \code{\link{copula-class}}.
}
\keyword{classes}
