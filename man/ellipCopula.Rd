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

  The  \code{\link{pCopula}()} methods for the normal- and t-copulas
  accept optional arguments to be passed to the underlying
  (numerical integration) algorithms from package \pkg{mvtnorm}'s
  \code{\link[mvtnorm]{pmvnorm}} and \code{\link[mvtnorm]{pmvt}},
  respectively, notably \code{algorithm}, see
  \code{\link[mvtnorm]{GenzBretz}}, or \code{abseps}
  which defaults to \code{0.001}.
## For smaller copula dimension 'd', alternatives are available and
## non-random, see ?GenzBretz from package 'mvtnorm' :



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
getSigma(t.cop)## the P ("Rho") matrix (with diagonal = 1)

## from the wrapper
norm.cop <- ellipCopula("normal", param = c(0.5, 0.6, 0.7),
                        dim = 3, dispstr = "un")
if(require("scatterplot3d") && dev.interactive(orNone=TRUE)) {
  ## 3d scatter plot of 1000 random observations
  scatterplot3d(rCopula(1000, norm.cop))
  scatterplot3d(rCopula(1000, t.cop))
}
set.seed(12)
uN <- rCopula(512, norm.cop)
set.seed(2); pN1 <- pCopula(uN, norm.cop)
set.seed(3); pN2 <- pCopula(uN, norm.cop)
stopifnot(all.equal(pN1, pN2, 1e-4))# see 5.711e-5
(Xtras <- copula:::doExtras())
if(Xtras) { ## a bit more accurately:
  set.seed(4); pN1. <- pCopula(uN, norm.cop, abseps = 1e-9)
  set.seed(5); pN2. <- pCopula(uN, norm.cop, abseps = 1e-9)
  stopifnot(all.equal(pN1., pN2., 1e-5))# see 3.397e-6
  ## but increasing the required precision (e.g., abseps=1e-15) does *NOT* help
}

## For smaller copula dimension 'd', alternatives are available and
## non-random, see ?GenzBretz from package 'mvtnorm' :
require("mvtnorm")# -> GenzBretz(), Miva(), and TVPACK() are available
## Note that Miwa() would become very slow for dimensions 5, 6, ..
set.seed(4); pN1.M <- pCopula(uN, norm.cop, algorithm = Miwa(steps = 512))
set.seed(5); pN2.M <- pCopula(uN, norm.cop, algorithm = Miwa(steps = 512))
stopifnot(all.equal(pN1.M, pN2.M, tol= 1e-15))# *no* randomness
set.seed(4); pN1.T <- pCopula(uN, norm.cop, algorithm = TVPACK(abseps = 1e-10))
set.seed(5); pN2.T <- pCopula(uN, norm.cop, algorithm = TVPACK(abseps = 1e-14))
stopifnot(all.equal(pN1.T, pN2.T, tol= 1e-15))# *no* randomness (but no effect of 'abseps')


## Versions with unspecified parameters:
tCopula()
allEQ <- function(u,v) all.equal(u, v, tolerance=0)
stopifnot(allEQ(ellipCopula("norm"), normalCopula()),
          allEQ(ellipCopula("t"), tCopula()))
tCopula(dim=3)
tCopula(dim=4, df.fixed=TRUE)
tCopula(dim=5, disp = "toep", df.fixed=TRUE)
normalCopula(dim=4, disp = "un")
}
\keyword{distribution}
\keyword{multivariate}
