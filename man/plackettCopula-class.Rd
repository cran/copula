\name{plackettCopula-class}
\docType{class}
\alias{plackettCopula-class}
\title{Class "plackettCopula" of Plackett Copulas}
\description{
  The Plackett copula class.
}
\section{Objects from the Class}{
  Objects can be created by \code{new("plackettCopula", ...)} but are
  typically produced by \code{\link{plackettCopula}(alpha)}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Numeric (scalar), the dimension of the copula.}
    \item{\code{exprdist}:}{a length two \code{\link{expression}} with
      expressions for the CDF and PDF of the copula.}% FIXME still used?
    \item{\code{parameters}:}{a number (numeric vector of length one) specifying the
      \emph{non negative} parameter.}
    \item{\code{param.names}:}{the \code{"\link{character}"} string \code{"alpha"}.}
    \item{\code{param.lowbnd}:}{the number \code{0}.}
    \item{\code{param.upbnd}:}{the number \code{Inf}.}
    \item{\code{fullname}:}{(deprecated; do not use!)}
  }
}
\section{Methods}{
  Typical copula methods work, see \code{"\linkS4class{plackettCopula}"}
  and use \code{methods(class = "plackettCopula")}.
}
\section{Extends}{
  Class \code{"plackettCopula"} extends class \code{"\linkS4class{copula}"} directly.
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Nelsen, R. B. (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{copula-class}}, \code{\link{plackettCopula}}.
}
\examples{
str(plackettCopula())

plackettCopula()@exprdist[["cdf"]]
methods(class = "plackettCopula")
contourplot2(plackettCopula(7), pCopula)
wireframe2(plackettCopula(5), dCopula, main= "plackettCopula(5)")
}
\keyword{classes}
