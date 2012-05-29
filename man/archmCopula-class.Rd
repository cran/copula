\name{archmCopula-class}
\docType{class}
\alias{archmCopula-class}
\alias{claytonCopula-class}
\alias{frankCopula-class}
\alias{gumbelCopula-class}
\alias{amhCopula-class}
\title{Class "archmCopula" }
\description{Archimedean copula class.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("archmCopula",
    ...)} or by function \code{\link{archmCopula}}. Implemented families are
  Clayton, Gumbel, Frank, and Ali-Mikhail-Haq.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"}:
      expressions of the cdf and pdf of the copula. These expressions
      are used in function 'pcopula' and 'dcopula'. }
    \item{\code{dimension}:}{Object of class \code{"numeric"}, dimension
    of the copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"},
      parameter values. }
    \item{\code{param.names}:}{Object of class \code{"character"},
      parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"},
      parameter lower bounds. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"},
      parameter upper bounds. }
    \item{\code{message}:}{Object of class \code{"character"}, family names
      of the copula. }
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
  Class \code{"archmCopula"} extends class \code{"\linkS4class{copula}"} directly.
  Class \code{"claytonCopula"}, \code{"frankCopula"},
  \code{"gumbelCopula"} and \code{"amhCopula"} extends class
  \code{"archmCopula"} directly.
}
\note{
  \code{"gumbelCopula"} is also of class \code{"\linkS4class{evCopula}"}.
}
\seealso{
  \code{\link{archmCopula}},
  \code{\link{copula-class}}.
}
\keyword{classes}
