\name{evCopula-class}
\docType{class}
\alias{evCopula-class}
\alias{galambosCopula-class}
\alias{huslerReissCopula-class}
\alias{tawnCopula-class}
\alias{tevCopula-class}
\title{Classes Representing Extreme-Value Copulas}
\description{
  Class \code{evCopula} is the virtual (mother) class of all
  extreme-value copulas.  There currently are five subclasses,
  \code{"galambosCopula"}, \code{"huslerReissCopula"},
  \code{"tawnCopula"}, \code{"tevCopula"}, and \code{"gumbelCopula"},
  the latter of which is also an Archimedean copula, see the page for
  class \code{"\linkS4class{archmCopula}"}.
}
\section{Objects from the Class}{
  \code{evCopula} is a virtual class: No objects may be created from it.
  Objects of class \code{"galambosCopula"} etc, can be created by calls
  of the form \code{new("galambosCopula", ...)}, but typically rather by
  \code{\link{galambosCopula}()}, etc, see there.
}
\section{Slots}{
  All slots are inherited from the mother class
  \code{"\linkS4class{copula}"}, see there.
}
\section{Methods}{
  \describe{
    \item{dCopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "galambosCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "huslerReissCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "huslerReissCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"evCopula"} extends class \code{"\linkS4class{copula}"}
  directly. Classes \code{"galambosCopula"}, \code{"huslerReissCopula"},
  \code{"tawnCopula"}, and \code{"tevCopula"} extend class
  \code{"evCopula"} directly.
}
\note{
  %The expressions of pdf are obtained by differentiating the cdf
  %expression using function \code{"deriv"}.
  Objects of class \code{"\linkS4class{gumbelCopula}"} are also of class
  \code{"\linkS4class{archmCopula}"}.
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{
  \code{\link{evCopula}},
  \code{\link{evTestC}},
  \code{\link{evTestK}},
  \code{\link{gofEVCopula}},
  \code{\link{copula-class}}.
}
\keyword{classes}
