\name{fgmCopula}
\alias{fgmCopula}

\title{ Construction of a fgmCopula Class Object}
\description{
  Constructs a multivariate multi-parameter Farlie-Gumbel-Morgenstern
  copula (see Nelsen 1999, p 87) class object with its corresponding
  parameters and dimension.
}
\usage{
fgmCopula(param, dim = 2)

}

\arguments{
  \item{param}{a numeric vector specifying the parameter values.}
  \item{dim}{the dimension of the copula.}
  \item{\dots}{ currently nothing. }
}
\value{
  A Farlie-Gumbel-Morgenstern copula object of class \code{"fgmCopula"}.
}
\note{
  The verification of the validity of the parameter values is of high
  complexity and may not work for high dimensional copulas. 

  The random number generation needs to be properly tested, especially
  for dimensions higher than 2.
}

\references{
  Nelsen (1999), \emph{An introduction to Copulas}, Springer, New York.
}

\author{Ivan Kojadinovic <ivan.kojadinovic@univ-nantes.fr> }

\seealso{ \code{\link{Copula}}, \code{\link{copula-class}}. }

\examples{
## a bivariate example
fgm.cop <- fgmCopula(1)
x <-rcopula(fgm.cop, 1000)
cor(x,method="kendall")
kendallsTau(fgm.cop)
cor(x,method="spearman")
spearmansRho(fgm.cop)
persp(fgm.cop, dcopula)
contour(fgm.cop, dcopula)

## a trivariate example with wrong parameter values
## fgm2.cop <- fgmCopula(c(1,1,1,1), dim=3)

## a trivariate example with satisfactory parameter values
fgm2.cop <- fgmCopula(c(.2,-1,-0.3,1), dim=3)
fgm2.cop
}
\keyword{ distribution }
\keyword{ multivariate }
