\name{Copula}
\alias{Copula}
\alias{dcopula}
\alias{pcopula}
\alias{rcopula}
\alias{dcopula,claytonCopula-method}
\alias{dcopula,frankCopula-method}
\alias{dcopula,gumbelCopula-method}
\alias{dcopula,normalCopula-method}
\alias{dcopula,tCopula-method}
\alias{pcopula,claytonCopula-method}
\alias{pcopula,frankCopula-method}
\alias{pcopula,gumbelCopula-method}
\alias{pcopula,normalCopula-method}
\alias{pcopula,tCopula-method}
\alias{rcopula,claytonCopula-method}
\alias{rcopula,frankCopula-method}
\alias{rcopula,gumbelCopula-method}
\alias{rcopula,normalCopula-method}
\alias{rcopula,tCopula-method}

\title{ The Copula Distribution }
\description{
  Density, distribution function, and random generation for a
  \code{"copula"} object.
}
\usage{
dcopula(copula, u)
pcopula(copula, u)
rcopula(copula, n)
}
\arguments{
  \item{copula}{ a \code{"copula"} object. }
  \item{u}{ a vector of the copula dimension or a matrix with number of
    rows being the copula dimension, giving the coordinates of the
    points where the density of distribution function need to be evaluated.
  }
  \item{n}{ number of observations to be generated. }
}
\details{
  The density function of an Archimedean copula is obtained by
  differentiating the distribution function symbolically using D.

  The distribution function of a t copula uses pmvt from package
  mvtnorm. The density function of a t copula uses the dmst from package
  sn.

  The random number generator for an Archiimedean copula uses the
  conditional approach for bivariate case and the Marshal-Olkin (1988)
  approach for dimension greater than 2.
}
\value{
  'dcopula' gives the density, 'pcopula' gives the distribution
  function, and 'rcopula' generates random variates.
}
\references{
  Joe (1997)
  Nelsen (1999)
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\seealso{ copulaObject, ellipCopula, archmCopula }
\examples{
norm.cop <- normalCopula(0.5)
norm.cop
x <- rcopula(norm.cop, 100)
plot(x)
dcopula(norm.cop, x)
pcopula(norm.cop, x)
persp(norm.cop, dcopula)
contour(norm.cop, pcopula)
## a 3-dimensional normal copula
u <- rcopula(normalCopula(0.5, dim = 3), 1000)
scatterplot3d(u)
## a 3-dimensional clayton copula
v <- rcopula(claytonCopula(2, dim = 3), 1000)
scatterplot3d(v)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
