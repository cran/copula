\name{plackettCopula}
\alias{plackettCopula}
\title{ Construction of Plackett Copula Class Object }
\description{
  Constructs an Plackett copula class object with its corresponding
  parameter.
}
\usage{
plackettCopula(param)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{ a numeric vector specifying the parameter values. }
}
\value{
  A Plackett copula object of class \code{"plackettCopula"}.
}
\author{ Jun Yan <jyan@stat.uconn.edu> }
\seealso{ \code{\link{ellipCopula}}, \code{\link{archmCopula}}  }
\examples{
plackett.cop <- plackettCopula(param=2)
tailIndex(plackett.cop)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
