\name{indepCopula-class}
\docType{class}
\alias{indepCopula-class}
\title{Class "indepCopula"}
\description{Independence copula class.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("indepCopula", ...)}
  or by function \code{\link{indepCopula}()}.  Such objects can be
  useful as special cases of parametric copulas, bypassing
  copula-specific computations such as distribution, density, and sampler.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Object of class \code{"numeric"}, dimension
      of the copula.}
    \item{\code{exprdist}:}{an \code{\link{expression}} of length two,
      for the \dQuote{formulas} of the cdf and pdf of the copula.}% used in ../R/dC-dc.R
  }
}
\section{Methods}{
  \describe{
    \item{A}{\code{signature(copula = "indepCopula")}: ... }
    \item{dCopula}{\code{signature(copula = "indepCopula")}: ... }
    \item{pCopula}{\code{signature(copula = "indepCopula")}: ... }
    \item{rCopula}{\code{signature(copula = "indepCopula")}: ... }
  }
}
\section{Extends}{
  Class \code{"indepCopula"} directly extends classes
  \code{"\linkS4class{dimCopula}"} and \code{"\linkS4class{parCopula}"}.
}
%\author{Ivan Kojadinovic, Jun Yan and Martin M}
\seealso{
  \code{\link{indepCopula}}; documentation for classes
  \code{\linkS4class{dimCopula}} and
  \code{\linkS4class{parCopula}}.
}
\examples{
getClass("indepCopula")
}
\keyword{classes}
