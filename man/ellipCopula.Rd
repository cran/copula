\name{ellipCopula}
\alias{ellipCopula}
\alias{normalCopula}
\alias{tCopula}
\title{ Construction of elliptical copula class object }
\description{
  Constructs an elliptical copula class object with its corresponding
  parameters and dimension.
}
\usage{
ellipCopula (family, param, dim = 2, dispstr = "ex", df = 4, ...)
normalCopula(param, dim = 2, dispstr = "ex")
    tCopula (param, dim = 2, dispstr = "ex", df = 4, df.fixed = FALSE)
}
\arguments{
  \item{family}{ a character string specifying the family of an
    elliptical copula. Implemented families are "normal" and "t".  }
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{ the dimension of the copula. }
  \item{dispstr}{ a character string specifying the type of the symmetric
    positive definite matrix characterizing the elliptical
    copula. Implemented structures are "ex" for exchangeable, "ar1" for
    AR(1), "toep" for Toeplitz, and "un" for unstructured. }
  \item{df}{ a integer value specifying the number of degrees of freedom of the
    multivariate t distribution used to construct the t copulas.}
  \item{df.fixed}{ TRUE means that the degrees of freedom will never be
    considered as a parameter to be estimated; FALSE means that
    \code{df} will be estimated if the object is passed as argument to \code{fitCopula}.}
  \item{\dots}{ currently nothing. }
}
\value{
  An elliptical copula object of class \code{"\linkS4class{normalCopula}"} or
  \code{"\linkS4class{tCopula}"}.
}
\note{
  \code{ellipCopula()} is a wrapper for \code{normalCopula()} and
  \code{tCopula()}.
}
\seealso{ \code{\link{archmCopula}}, \code{\link{fitCopula}}.
}
\examples{
norm.cop <- normalCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un")
t.cop <- tCopula(c(0.5, 0.3), dim = 3, dispstr = "toep",
                 df = 2, df.fixed = TRUE)
## from the wrapper
norm.cop <- ellipCopula("normal", param = c(0.5, 0.6, 0.7),
                        dim = 3, dispstr = "un")
if(require("scatterplot3d")) {
  ## 3d scatter plot of 1000 random observations
  scatterplot3d(rcopula(norm.cop, 1000))
  scatterplot3d(rcopula(t.cop, 1000))
}
}
\keyword{distribution}
\keyword{multivariate}
