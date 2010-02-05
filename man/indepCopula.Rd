\name{indepCopula}
\alias{indepCopula}
\title{Construction of independence copula class objects}
\description{
  Constructs an independence copula class object with its corresponding
  dimension.
}
\usage{
indepCopula(dim = 2)
}
\arguments{
  \item{dim}{ the dimension of the copula. }
}
\value{
  An independence copula object of class \code{"indepCopula"}.
}
\seealso{ \code{\link{archmCopula}}, \code{\link{ellipCopula}}, \code{\link{evCopula}}. }
\examples{
indep.cop <- indepCopula(3)
x <- rcopula(indep.cop, 10)
dcopula(indep.cop, x)
persp(indepCopula(), pcopula)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
