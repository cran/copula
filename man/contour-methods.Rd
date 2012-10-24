\name{contour-methods}
\docType{methods}
\alias{contour-methods}
\alias{contour,indepCopula-method}
\alias{contour,copula-method}
\alias{contour,mvdc-method}
\title{Methods for Function `contour' in Package `copula'}
\description{
 Methods for function \code{contour} in package \pkg{copula}
}
\section{Methods}{
\describe{
  \item{x = "copula"}{ Contour plot for a \code{"\linkS4class{copula}"} or
    \code{"\linkS4class{indepCopula}"} object. }

  \item{x = "mvdc"}{ Contour plot for a \code{"mvdc"} object. }
}}
\details{
  When \code{x} is of class \code{"\linkS4class{copula}"}, the following
  arguments can be supplied:
  \describe{
    \item{fun}{the function to be plotted, i.e.,\code{\link{dCopula}} or
      \code{\link{pCopula}}.}
    \item{n (= 51)}{the number of points to do the plotting.}
    \item{arguments for \code{"contour"}}{theta=-30, phi=30, expand=0.618.}
  }
  When \code{x} is of class \code{"mvdc"}, the following arguments
  replace the effect of \code{n = 51}:
  \describe{
    \item{xlim}{the range of the x variable.}
    \item{ylim}{the range of the y variable.}
    \item{nx}{the number of points for x to expand.}
    \item{ny}{the number of points for y to expand.}
  }
}
\seealso{
  The \code{\link{persp-methods}} for \dQuote{perspective} aka
  \dQuote{3D} plots.
}
\examples{
contour(frankCopula(-0.8), dCopula)
contour(frankCopula(-1.2), pCopula)
contour(claytonCopula(2), pCopula)
x <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd =1), list(mean = 1)))
contour(x, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
contour(x, pMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
