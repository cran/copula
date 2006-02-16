\name{ellipCopula}
\alias{ellipCopula}
\alias{normalCopula}
\alias{tCopula}
\title{ Construction of Elliptical Copula Class Object }
\description{
  Constructs an elliptical copula class object with its corresponding
  parameter and dimension.
}
\usage{
ellipCopula(family, param, dim = 2, dispstr = "ex", df = 5, ...)
normalCopula(param, dim = 2, dispstr = "ex")
tCopula(param, dim = 2, dispstr = "ex", df = 5)
}
\arguments{
  \item{family}{ a character string specifying the family of an
    elliptical copula. Implemented families are "normal" and "t".  }
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{ the dimension of the copula. }
  \item{dispstr}{ a character string specifying the type of the symmetric
    positive definite matrix characterizing the elliptical
    copula. Implemented structures are "ex" for exchangeable, "ar1" for
    AR(1), "toep" for toeplitz, , and "un" for unstructured. For normal
    copula, this defines the structure of the correlation matrix. }
  \item{df}{ a numerical value specifying the degree of freedom for the
    multivariate t distribution used to construct the t copulas.}
  \item{\dots}{ currently nothing. }
}
\value{
  An elliptical copula object of class \code{"normalCopula"} or
  \code{"tCopula"}.
}
\note{
  \code{"ellipCopula"} is a wrapper for \code{"normalCopula"} and
  \code{"tCopula"}.
}
\author{ Jun Yan <jyan@stat.uiowa.edu> }
\seealso{ \code{\link{archmCopula}} }
\examples{
norm.cop <- normalCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un")
t.cop <- tCopula(c(0.5, 0.3), dim = 3, dispstr = "toep", df = 2)
## from the wrapper
norm.cop <- ellipCopula("normal", param = c(0.5, 0.6, 0.7),
                        dim = 3, dispstr = "un")
## 3d scatter plot of 1000 random observations
##scatterplot3d(rcopula(norm.cop, 1000))
##scatterplot3d(rcopula(t.cop, 1000))
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
