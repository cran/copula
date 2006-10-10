\name{archmCopula}
\alias{archmCopula}
\alias{claytonCopula}
\alias{frankCopula}
\alias{gumbelCopula}
\title{ Construction of Archimedean Copula Class Object }
\description{
  Constructs an Archimedean copula class object with its corresponding
  parameter and dimension.
}
\usage{
archmCopula(family, param, dim = 2, ...)
claytonCopula(param, dim = 2)
frankCopula(param, dim = 2)
gumbelCopula(param, dim = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{ a character string specifying the family of an
    Archimedean copula. Implemented families are "clayton", "frank", and
    "gumbel".
  }
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{ the dimension of the copula. }
  \item{\dots}{ currently nothing. }
}
\value{
  An elliptical copula object of class \code{"claytonCopula"},
  \code{"frankCopula"}, or \code{"gumbelCopula"}.
}
\note{
  \code{"archmCopula"} is a wrapper for \code{"claytonCopula"},
  \code{"frankCopula"}, and \code{"gumbelCopula"}.

  For dim > 6, the expression of pdf is not available due to intensive
  compting involved in symbolically differentiating the cdf. Therefore,
  for dim > 6, likelihood estimation cannot be done yet.
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\seealso{ \code{\link{ellipCopula}} }
\examples{
clayton.cop <- claytonCopula(2, dim = 3)
##scatterplot3d(rcopula(clayton.cop, 1000))
frank.cop <- frankCopula(3)
persp(frank.cop, dcopula)
gumbel.cop <- archmCopula("gumbel", 5)
contour(gumbel.cop, dcopula)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
