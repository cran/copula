\name{fgmCopula-class}
\docType{class}
\alias{fgmCopula-class}
\title{Class "fgmCopula" - Multivariate Multiparameter Farlie-Gumbel-Morgenstern Copulas}
\description{
  The class of multivariate multiparameter Farlie-Gumbel-Morgenstern
  copulas are typically created via \code{\link{fgmCopula}(..)}.
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{fgmCopula}(..)}, or more
  low-level by (careful) calls to \code{new("fgmCopula", ..)}.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"},
      expressions for the cdf and pdf of the copula. These expressions
      are used in function \code{pCopula()} and \code{dCopula()}.}
    \item{\code{subsets.char}:}{Object of class \code{"character"},
      containing the subsets of integers used for naming the parameters. }
    \item{\code{dimension}:}{Object of class \code{"numeric"}, the
      dimension of the copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"},
      parameter values. }
    \item{\code{param.names}:}{Object of class \code{"character"},
      parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"},
      parameter lower bound. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"},
      parameter upper bound. }
    \item{\code{fullname}:}{Object of class \code{"character"}, family names
      of the copula (deprecated).}
  }
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "fgmCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"fgmCopula"} extends class \code{"\linkS4class{copula}"} directly.
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Nelsen, R. B. (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\note{
  The verification of the validity of the parameter values is of high
  complexity and may not work for high dimensional copulas.

  The random number generation needs to be properly tested, especially
  for dimensions higher than 2.
}
\seealso{
  \code{\link{copula-class}}, \code{\link{fgmCopula-class}}; to create such
  objects, use \code{\link{fgmCopula}()}; see there, also for examples.
}
\keyword{classes}
