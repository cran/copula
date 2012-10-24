\name{mvdc-class}
\title{Class "mvdc"}
\docType{class}
\alias{mvdc-class}
\alias{dim,mvdc-method}
\alias{show,mvdc-method}
\description{Class representing multivariate distributions constructed
  using Sklar's theorem. }
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mvdc}()}, or
  can be created by calls of the form \code{new("mvdc", ...)}.
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
    \item{dim}{\code{signature(x = "mvdc")}: the dimension of the
      distribution; this is the same as \code{dim(x@copula)}.}
    \item{persp}{\code{signature(x = "mvdc")}: ... }
    \item{show}{\code{signature(object = "mvdc")}: quite compactly display
      the content of the "mvdc" \code{object}.}
  }
}
\seealso{
  \code{\link{mvdc}}, % ./Mvdc.Rd
  also for examples.
}
\keyword{classes}
