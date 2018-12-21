\name{indepCopula}
\alias{indepCopula}
\title{Construction of Independence Copula Objects}
\description{
  Constructs an independence copula with its corresponding dimension.
}
\usage{
indepCopula(dim = 2)
}
\arguments{
  \item{dim}{the dimension of the copula.}
}
\value{
  An independence copula object of class \code{"\linkS4class{indepCopula}"}.
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{Mathematically, the independence copula is also a special
  (boundary) case of e.g., classes \code{"\linkS4class{archmCopula}"},
  \code{"\linkS4class{ellipCopula}"}, and \code{"\linkS4class{evCopula}"}.
}
\examples{
indep.cop <- indepCopula(3)
x <- rCopula(10, indep.cop)
dCopula(x, indep.cop)
persp(indepCopula(), pCopula)
}
\keyword{ distribution }
\keyword{ multivariate }
