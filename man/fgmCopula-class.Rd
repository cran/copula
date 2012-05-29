\name{fgmCopula-class}
\docType{class}
\alias{fgmCopula-class}

\title{Class "fgmCopula" }
\description{Multivariate multiparameter Farlie-Gumbel-Morgenstern
  copula.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fgmCopula",
    ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"},
      expressions for the cdf and pdf of the copula. These expressions
      are used in function \code{"pcopula"} and \code{"dcopula"}. }
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
    \item{\code{message}:}{Object of class \code{"character"}, family names
      of the copula. }
  }
}
\section{Methods}{
  \describe{
    \item{dcopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "fgmCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"fgmCopula"} extends class \code{"\linkS4class{copula}"} directly.
}

\references{
  R.B. Nelsen (2006), \emph{An introduction to Copulas}, Springer, New York.
}

\note{
  The verification of the validity of the parameter values is of high
  complexity and may not work for high dimensional copulas.

  The random number generation needs to be properly tested, especially
  for dimensions higher than 2.
}
\seealso{
  \code{\link{copula-class}}, \code{\link{fgmCopula-class}}.
}
\keyword{classes}
