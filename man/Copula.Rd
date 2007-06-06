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
\alias{dcopula,galambosCopula-method}
\alias{dcopula,amhCopula-method}
\alias{dcopula,huslerReissCopula-method}
\alias{dcopula,plackettCopula-method}
\alias{dcopula,fgmCopula-method} 

\alias{pcopula,claytonCopula-method}
\alias{pcopula,frankCopula-method}
\alias{pcopula,gumbelCopula-method}
\alias{pcopula,normalCopula-method}
\alias{pcopula,tCopula-method}
\alias{pcopula,galambosCopula-method}
\alias{pcopula,amhCopula-method}
\alias{pcopula,huslerReissCopula-method}
\alias{pcopula,plackettCopula-method}
\alias{pcopula,fgmCopula-method} 

\alias{rcopula,claytonCopula-method}
\alias{rcopula,frankCopula-method}
\alias{rcopula,gumbelCopula-method}
\alias{rcopula,normalCopula-method}
\alias{rcopula,tCopula-method}
\alias{rcopula,galambosCopula-method}
\alias{rcopula,amhCopula-method}
\alias{rcopula,huslerReissCopula-method}
\alias{rcopula,plackettCopula-method}
\alias{rcopula,fgmCopula-method} 

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
  Joe (1997), \emph{Multivariate Models and Dependence Concepts},
  Chapman and Hall, London.
  
  Nelsen (1999), \emph{An introduction to Copulas}, Springer, New York.
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\seealso{ \code{\link{copula-class}}, \code{\link{ellipCopula}}, \code{\link{archmCopula}}, \code{\link{fgmCopula}}. }
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
## scatterplot3d(u)
## a 3-dimensional clayton copula
v <- rcopula(claytonCopula(2, dim = 3), 1000)
## scatterplot3d(v)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
