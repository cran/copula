\name{copula-class}
\docType{class}
\alias{copula-class}
\alias{plackettCopula-class}
\title{Class "copula" }
\description{ A class of multivariate distribution with uniform margins. }
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
    \item{\code{message}:}{Object of class \code{"character"}, families
      of the copula. }
  }
}
\author{ Jun  Yan <jyan@stat.uiowa.edu> }
\note{
  The 'copula' class is extended by the 'evCopula', 'archmCopula' and
  'ellipCopula' classes. Objects of implemented copulas can be created
  from functions 'evCopula', 'archmCopula' and 'ellipcopula'.

  \code{"plackettCopula"} is a special type of copula which does not
  belong to either one of the three classes above.
}
\section{Warning }{ This implementation is still at the experimental
  stage and is subject to change during the development.
} 
\seealso{
  \code{\link{archmCopula-class}},
  \code{\link{ellipCopula-class}}
  \code{\link{evCopula-class}}
}
\keyword{classes}
