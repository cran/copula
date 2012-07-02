\name{mvdc-class}
\docType{class}
\alias{mvdc-class}
\title{Class "mvdc"}
\description{ Class representing multivariate distributions constructed
  using Sklar's theorem. }
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("mvdc", ...)}
  or by function \code{mvdc}.
}
\section{Slots}{
  \describe{
    \item{\code{copula}:}{Object of class \code{"\linkS4class{copula}"},
      specifying the copula.}
    \item{\code{margins}:}{Object of class \code{"character"},
      specifying the marginal distributions. }
    \item{\code{paramMargins}:}{Object of class \code{"list"}, whose
      each component is a list of named components, giving the parameter
      values of the marginal distributions. See \code{\link{mvdc}}. }
    \item{\code{marginsIdentical}:}{Object of class \code{"logical"},
      that, if TRUE, restricts the marginal distributions to be
      identical, default is \code{FALSE}. }
  }
}
\section{Methods}{
  \describe{
    \item{contour}{\code{signature(x = "mvdc")}: ... }
    \item{persp}{\code{signature(x = "mvdc")}: ... }
  }
}
\seealso{
  \code{\link{mvdc}}.
}
\keyword{classes}
