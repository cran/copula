\name{archmCopula-class}
\docType{class}
\alias{archmCopula-class}
\alias{claytonCopula-class}
\alias{frankCopula-class}
\alias{gumbelCopula-class}
\alias{amhCopula-class}
\title{Class "archmCopula" }
\description{ Archimedean copula }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("archmCopula",
    ...)} or by function 'archmCopula'.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"},
      expressions for the cdf and pdf of the copula. These expressions
      are used in function 'pcopula' and 'dcopula'. }
    \item{\code{dimension}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{parameters}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.names}:}{Object of class \code{"character"} ~~ }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"} ~~ }
    \item{\code{message}:}{Object of class \code{"character"} ~~ }
  }
}
\section{Methods}{
  \describe{
    \item{dcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "amhCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"archmCopula"} extends class \code{"copula"} directly.
  Class \code{"claytonCopula"}, \code{"frankCopula"},
  \code{"gumbelCopula"} and \code{"amhCopula"} extends class
  \code{"archmCopula"} directly.
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\note{
  The expressions of pdf are obtained by differentiating the cdf
  expression using function 'D'.

  \code{"gumbelCopula"} is also of class \code{"evCopula"}.
}
\seealso{
  \code{\link{archmCopula}}, 
  \code{\link{copula-class}}.
}
\keyword{classes}
