\name{fgmCopula}
\alias{fgmCopula}

\title{ Construction of a fgmCopula class object}
\description{
  Constructs a multivariate multiparameter Farlie-Gumbel-Morgenstern
  copula class object with its corresponding parameters and dimension.
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
  R.B. Nelsen (2006), \emph{An introduction to Copulas}, Springer, New York.
}

\seealso{
  \code{\link{Copula}}, \code{\link{copula-class}},
  \code{\link{fitCopula}}.
}
\examples{
## a bivariate example
fgm.cop <- fgmCopula(1)
x <- rcopula(fgm.cop, 1000)
cor(x, method = "kendall")
kendallsTau(fgm.cop)
cor(x, method = "spearman")
spearmansRho(fgm.cop)
persp(fgm.cop, dcopula)
contour(fgm.cop, dcopula)

## a trivariate example with wrong parameter values
## fgm2.cop <- fgmCopula(c(1,1,1,1), dim = 3)

## a trivariate example with satisfactory parameter values
fgm2.cop <- fgmCopula(c(.2,-.2,-.4,.6), dim = 3)
fgm2.cop
}
\keyword{ distribution }
\keyword{ multivariate }
