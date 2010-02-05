\name{evCopula-class}
\docType{class}
\alias{evCopula-class}
\alias{galambosCopula-class}
\alias{huslerReissCopula-class}
\alias{tawnCopula-class}
\alias{tevCopula-class}

\title{Class "evCopula" }
\description{Class representing extreme-value copulas.}
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
  Class \code{"evCopula"} extends class \code{"copula"}
  directly. Classes \code{"galambosCopula"}, \code{"huslerReissCopula"},
  \code{"tawnCopula"}, and \code{"tevCopula"} extend class \code{"evCopula"} directly.
}
\note{
  %The expressions of pdf are obtained by differentiating the cdf
  %expression using function \code{"deriv"}.
  Objects of class \code{"gumbelCopula"} are also of class \code{"archmCopula"}.
}
\seealso{
  \code{\link{evCopula}},
  \code{\link{evTestC}},
  \code{\link{evTestK}},
  \code{\link{gofEVCopula}},
  \code{\link{copula-class}}.
}
\keyword{classes}
