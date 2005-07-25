\name{contour-methods}
\docType{methods}
\alias{contour-methods}
\alias{contour,copula-method}
\alias{contour,mvdc-method}
\title{ Methods for Function contour in Package `copula' }
\description{
 Methods for function \code{contour} in package \pkg{copula} 
}
\section{Methods}{
\describe{
\item{x = "copula"}{ Contour plot for a \code{"copula"} object. }

\item{x = "mvdc"}{ Contour plot for a \code{"mvdc"} object. }
}}
\details{
  When x is of class \code{"copula"}, these arguments can be supplied:
  fun: the function to be plotted, \code{"dcopula"} or \code{"pcopula"}.
  n = 51: the number of points to do the plotting.
  theta = -30, phi = 30, expand = 0.618: arguments for \code{"contour"}

  when x is of class \code{"mvdc"}, two more arguments are expected to
  replace the effect of n = 51:
  xis: the values of the x variable.
  yis: the values of the y variable.
}
\examples{
contour(frankCopula(-0.8), dcopula)
contour(claytonCopula(2), pcopula)
x <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd =1), list(mean = 1)))
contour(x, dmvdc, xis = seq(-2, 2, len = 51), yis = seq(-1, 3, len = 51))
contour(x, pmvdc, xis = seq(-2, 2, len = 51), yis = seq(-1, 3, len = 51))
}
\keyword{methods}
\keyword{hplot}
