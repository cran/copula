\name{archmCopula-class}
\docType{class}
\alias{archmCopula-class}
\alias{claytonCopula-class}
\alias{frankCopula-class}
\alias{gumbelCopula-class}
\alias{amhCopula-class}
\alias{joeCopula-class}
\title{Class "archmCopula"}
\description{Archimedean copula class.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("archmCopula",
    ...)} or by function \code{\link{archmCopula}}.  Implemented
  families are Clayton, Gumbel, Frank, Joe, and Ali-Mikhail-Haq.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"}:
      expressions of the cdf and pdf of the copula. These expressions
      are used in function \code{\link{pCopula}} and \code{\link{dCopula}}. }
    \item{\code{dimension, parameters}, etc:}{all inherited from the
      super class \code{\linkS4class{copula}}.}
  }
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "claytonCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "frankCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "gumbelCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "amhCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "joeCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "joeCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "joeCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"archmCopula"} extends class \code{"\linkS4class{copula}"}
  directly.  Class \code{"claytonCopula"}, \code{"frankCopula"},
  \code{"gumbelCopula"}, \code{"amhCopula"} and \code{"joeCopula"}
  extends class \code{"archmCopula"} directly.
}
\note{
  \code{"gumbelCopula"} is also of class \code{"\linkS4class{evCopula}"}.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\seealso{
  \code{\link{archmCopula}}, for constructing such copula objects;
  \code{\link{copula-class}}.
}
\keyword{classes}
