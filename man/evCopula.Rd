\name{evCopula}
\alias{evCopula}
\alias{galambosCopula}
\alias{huslerReissCopula}
%\alias{gumbelCopula}
\title{ Construction of Extreme Value Copula Class Object }
\description{
  Constructs an extreme value copula class object with its corresponding
  parameter.
}
\usage{
evCopula(family, param, dim = 2, ...)
galambosCopula(param)
huslerReissCopula(param)
%gumbelCopula(param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{family}{ a character string specifying the family of an
    extreme value copula. Implemented families are "galambos" and
    "gumbel".
  }
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{ the dimension of the copula. }
  \item{\dots}{ currently nothing. }
}
\value{
  An extreme value copula object of class \code{"galambosCopula"}
  \code{"gumbelCopula"}, or \code{"huslerReissCopula"}.
}
\note{
  Gumbel copula is both Archimedean copula and extreme value copula.
}
\author{ Jun Yan <jyan@stat.uconn.edu> }
\seealso{ \code{\link{ellipCopula}}, \code{\link{archmCopula}}  }
\examples{
gumbel.cop <- evCopula("gumbel", param=2)
contour(gumbel.cop, dcopula)
galambos.cop <- galambosCopula(2)
contour(galambos.cop, dcopula)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
