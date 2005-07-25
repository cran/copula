\name{Mvdc}
\alias{Mvdc}
\alias{mvdc}
\alias{dmvdc}
\alias{pmvdc}
\alias{rmvdc}
\title{ Multivariate Distribution via Copula }
\description{
  Density, distribution function, and random generator for a
  multivariate distribution via copula.
}
\usage{
mvdc(copula, margins, paramMargins)
dmvdc(mvdc, x)
pmvdc(mvdc, x)
rmvdc(mvdc, n)
}
\arguments{
  \item{copula}{ an object of copula. }
  \item{margins}{ a character vector specifying all the marginal
    distributions. See details below.}
  \item{paramMargins}{ a list of list with names components, giving the
    parameter values of the marginal distributions. See details below. }
  \item{mvdc}{ a \code{mvdc} object. }
  \item{x}{ a vector of the copula dimension or a matrix with number of
    rows being the copula dimension, giving the coordinates of the
    points where the density of distribution function need to be evaluated.
  }
  \item{n}{ number of observations to be generated. }
}
\details{
  The characters in argument \code{margins} are used to construct
  function names of density, distribution, and quantile functions. For
  example, "norm" can be used to specify marginal distribution, because
  "dnorm", "pnorm", and "qnorm" are all available.

  Each component list in argument \code{paramMargins} is a list with
  named component which are used to specify the parameters of the
  marginal distributions. For example,
  \code{paramMargins = list(list(mean = 0, sd = 2), list(rate = 2))}
  can be used to specify that the first margin is normal with mean 0 and
  sd 2, and the second margin is exponential with rate 2.
}
\value{
  'mvdc' constructs an object of class \code{"mvdc"}.
  'dmvdc' gives the density, 'pmvdc' gives the distribution
  function, and 'rmvdc' generates random variates.
}
\author{  Jun Yan <jyan@stat.uiowa.edu> }
\seealso{
  \code{\link{ellipCopula}},
  \code{\link{archmCopula}},
  \code{\link{mvdc-class}},
  \code{\link{copula-class}}
}
\examples{
## construct a bivariate distribution whose marginals
## are normal and exponential respectively, coupled
## together via a normal copula
x <- mvdc(normalCopula(0.75), c("norm", "exp"),
          list(list(mean = 0, sd =2), list(rate = 2)))
x.samp <- rmvdc(x, 100)
dmvdc(x, x.samp)
pmvdc(x, x.samp)
persp(x, dmvdc, xis = seq(-4, 4, len = 51), yis = seq(0, 0.5, len = 51))
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
