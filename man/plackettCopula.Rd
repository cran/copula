\name{plackettCopula}
\alias{plackettCopula}
\title{Construction of a Plackett copula class object}
\description{
  Constructs a Plackett copula class object with its corresponding
  parameter.
}
\usage{
plackettCopula(param)
}
\arguments{
  \item{param}{ a numeric vector specifying the parameter values. }
}
\value{
  A Plackett copula object of class \code{"plackettCopula"}.
}
\references{
  R.L. Plackett (1965). A Class of Bivariate
  Distributions. \emph{Journal of the American Statistical Association},
  60:516--522.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{archmCopula}}.
}
\examples{
plackett.cop <- plackettCopula(param=2)
tailIndex(plackett.cop)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
