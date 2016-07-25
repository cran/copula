\name{indepCopula}
\alias{indepCopula}
\title{Construction of Independence Copula Class Objects}
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
  An independence copula object of class \code{"\linkS4class{indepCopula}"}.
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{\code{\link{archmCopula}}, \code{\link{ellipCopula}},
  \code{\link{evCopula}}.}
\examples{
indep.cop <- indepCopula(3)
x <- rCopula(10, indep.cop)
dCopula(x, indep.cop)
persp(indepCopula(), pCopula)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }
