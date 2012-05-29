\name{persp-methods}
\docType{methods}
\alias{persp-methods}
\alias{persp,indepCopula-method}
\alias{persp,copula-method}
\alias{persp,mvdc-method}
\title{ Methods for function `persp' in Package `copula' }
\description{
 Methods for function \code{persp} in package \pkg{copula}
}
\section{Methods}{
\describe{
  \item{x = "copula"}{perspective plot for a \code{"\linkS4class{copula}"} or
    \code{"\linkS4class{indepCopula}"} object. }
  \item{x = "mvdc"}{perspective plot for a \code{"mvdc"} object. }
}}
\details{
  When \code{x} is of class \code{"\linkS4class{copula}"}, the following arguments can be supplied:
  \describe{
    \item{fun}{the function to be plotted, i.e., \code{"dcopula"} or \code{"pcopula"}.}
    \item{n (= 51)}{the number of points to do the plotting.}
    \item{arguments for \code{"persp"}}{theta = -30, phi = 30, expand = 0.618.}
  }

  When \code{x} is of class \code{"mvdc"}, the following arguments
  replace the effect of \code{n = 51}:
  \describe{
    \item{xlim}{the range of the x variable.}
    \item{ylim}{the range of the y variable.}
  }
}
\examples{
persp(frankCopula(-0.8), dcopula)
persp(claytonCopula(2), pcopula)
x <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd =1), list(mean = 1)))
persp(x, dmvdc, xlim=c(-2, 2), ylim=c(-1, 3))
persp(x, pmvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
