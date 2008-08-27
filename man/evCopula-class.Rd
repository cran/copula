\name{evCopula-class}
\docType{class}
\alias{evCopula-class}
\alias{galambosCopula-class}
\alias{huslerReissCopula-class}
\title{Class "evCopula" }
\description{Class representing extreme value copulas.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("evCopula",
    ...)} or by function \code{"evCopula"}.
}
\section{Slots}{
  \describe{
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
    \item{dcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"evCopula"} extends class \code{"copula"} directly. Both
  class \code{"galambosCopula"} and class \code{"huslerReissCopula"}
  extend class \code{"evCopula"} directly.
}
\note{
  The expressions of pdf are obtained by differentiating the cdf
  expression using function \code{"deriv"}.

  \code{"gumbelCopula"} is also of class \code{"archmCopula"}.
}
\seealso{
  \code{\link{evCopula}}, 
  \code{\link{copula-class}}.
}
\keyword{classes}
