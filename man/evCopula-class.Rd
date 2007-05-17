\name{evCopula-class}
\docType{class}
\alias{evCopula-class}
\alias{galambosCopula-class}
\alias{huslerReissCopula-class}
\title{Class "evCopula" }
\description{ Extreme value copula }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("evCopula",
    ...)} or by function 'evCopula'.
}
\section{Slots}{
  \describe{
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
    \item{dcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "huslerReissCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"evCopula"} extends class \code{"copula"} directly.
  Class \code{"galambosCopula"} and \code{"huslerReissCopula"}
  extends class \code{"evCopula"} directly.
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\note{
  The expressions of pdf are obtained by differentiating the cdf
  expression using function 'deriv'.

  \code{"gumbelCopula"} is also of class \code{"archmCopula"}.
}
\seealso{
  \code{\link{evCopula}}, 
  \code{\link{copula-class}}.
}
\keyword{classes}
