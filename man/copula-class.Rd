\name{copula-class}
\docType{class}
\alias{copula-class}
\alias{plackettCopula-class}
\title{Class "copula" }
\description{ A class representing multivariate distributions with uniform margins. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("copula", ...)}.
}
\section{Slots}{
  \describe{
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
\note{
  The \code{"copula"} class is extended by the \code{"evCopula"},
  \code{"archmCopula"}, and \code{"ellipCopula"} classes. Instances of
  implemented copulas can be created from functions \code{"evCopula"},
  \code{"archmCopula"} and \code{"ellipcopula"}.
  
  \code{"plackettCopula"} and \code{"fgmCopula"} are special types of
  copulas which do not belong to either one of the three classes above.
}
\section{Warning }{
  This implementation is still at the experimental
  stage and is subject to change during the development.
} 
\seealso{
  \code{\link{archmCopula-class}},
  \code{\link{ellipCopula-class}},
  \code{\link{evCopula-class}},
  \code{\link{fgmCopula-class}}.
}
\keyword{classes}
