\name{ellipCopula}
\alias{ellipCopula}
\alias{normalCopula}
\alias{tCopula}
\title{Construction of Elliptical Copula Class Object}
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
  \item{param}{a numeric vector specifying the parameter values.  The
    \code{getRho()} method accesses this vector, whereas
    \code{\link{p2P}()} and \code{\link{getSigma}()} provide
    the corresponding \dQuote{Rho} matrix, see below.}
  \item{dim}{ the dimension of the copula. }
  \item{dispstr}{ a character string specifying the type of the symmetric
    positive definite matrix characterizing the elliptical
    copula. Implemented structures are "ex" for exchangeable, "ar1" for
    AR(1), "toep" for Toeplitz, and "un" for unstructured. }
  \item{df}{ a integer value specifying the number of degrees of freedom
    of the multivariate t distribution used to construct the t copulas.}
  \item{df.fixed}{logical specifying if the degrees of freedom \code{df} will be
    considered as a parameter (to be estimated) or not.  The default,
    \code{FALSE}, means that \code{df} is to be estimated if the object is
    passed as argument to \code{fitCopula}.}
  \item{\dots}{ currently nothing. }
}
\value{
  An elliptical copula object of class \code{"\linkS4class{normalCopula}"}
  or \code{"\linkS4class{tCopula}"}.
}
\note{
  \code{ellipCopula()} is a wrapper for \code{normalCopula()} and
  \code{tCopula()}.
}
\seealso{
  \code{\link{p2P}()}, and \code{\link{getSigma}()} for construction and
  extraction of the \dQuote{Rho} (\eqn{P}) or \eqn{Sigma} matrix of
  (generalized)% for t
  correlations.

  \code{\link{archmCopula}}, \code{\link{fitCopula}}.
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
  scatterplot3d(rCopula(1000, norm.cop))
  scatterplot3d(rCopula(1000, t.cop))
}

## Versions with unspecified parameters:
tCopula()
allEQ <- function(u,v) all.equal(u, v, tol=0)
stopifnot(allEQ(ellipCopula("norm"), normalCopula()),
          allEQ(ellipCopula("t"), tCopula()))
tCopula(dim=3)
tCopula(dim=4, df.fixed=TRUE)
tCopula(dim=5, disp = "toep", df.fixed=TRUE)
normalCopula(dim=4, disp = "un")
}
\keyword{distribution}
\keyword{multivariate}
