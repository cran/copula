\name{mvdc-class}
\docType{class}
\alias{mvdc-class}
\title{Class "mvdc" }
\description{ A class of multivariate distribution via copula. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("mvdc", ...)}
  or by function \code{mvdc}.
}
\section{Slots}{
  \describe{
    \item{\code{copula}:}{Object of class \code{"copula"}, specifying
      the copula. }
    \item{\code{margins}:}{Object of class \code{"character"},
      specifying the marginal distributions. }
    \item{\code{paramMargins}:}{Object of class \code{"list"}, a list of
    list, with eaching list giving the names parameter values of the margins.}
  }
}
\section{Methods}{
  \describe{
    \item{contour}{\code{signature(x = "mvdc")}: ... }
    \item{persp}{\code{signature(x = "mvdc")}: ... }
  }
}
\author{ Jun Yan <jyan@stat.uconn.edu> }
\seealso{
  \code{\link{mvdc}}
}
\keyword{classes}
