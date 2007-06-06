\name{fgmCopula-class}
\docType{class}
\alias{fgmCopula-class}

\title{Class "fgmCopula" }
\description{Multivariate multi-parameter Farlie-Gumbel-Morgenstern
  copula as defined in Nelsen (1999, p 87).}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("fgmCopula",
    ...)}.
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
    \item{dcopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "fgmCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "fgmCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"fgmCopula"} extends class \code{"copula"} directly.
}

\references{
  Nelsen (1999), \emph{An introduction to Copulas}, Springer, New York.
}

\author{Ivan Kojadinovic <ivan.kojadinovic@univ-nantes.fr> }

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
