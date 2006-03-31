\name{persp-methods}
\docType{methods}
\alias{persp-methods}
\alias{persp,copula-method}
\alias{persp,mvdc-method}
\title{ Methods for Function persp in Package `copula' }
\description{
 Methods for function \code{persp} in package \pkg{copula} 
}
\section{Methods}{
\describe{
  \item{x = "copula"}{ Perspective plot for a \code{"copula"} object. }
  \item{x = "mvdc"}{ Perspective plot for a \code{"mvdc"} object. }
}}
\details{
  When x is of class \code{"copula"}, these arguments can be supplied:
  fun: the function to be plotted, \code{"dcopula"} or \code{"pcopula"}.
  n = 51: the number of points to do the plotting.
  theta = -30, phi = 30, expand = 0.618: arguments for \code{"persp"}

  when x is of class \code{"mvdc"}, two more arguments are expected to
  replace the effect of n = 51:
  xlim: the range of the x variable.
  ylim: the range of the y variable.
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
