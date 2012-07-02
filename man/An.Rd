\name{An}
\alias{An.biv}
\alias{An}
\alias{Anfun}%deprecated
\title{Nonparametric rank-based estimators of Pickands dependence function}
\description{
  Bivariate and multivariate versions of the nonparametric rank-based
  estimators of the Pickands dependence function \eqn{A}, studied in
  Genest and Segers (2009) and Gudendorf and Segers (2011).
}
\usage{
An.biv(x, w, estimator = c("CFG", "Pickands"), corrected = TRUE)
An(x, w)
}
\arguments{
  \item{x}{ a data matrix that will be transformed to
    pseudo-observations. If \code{An.biv} is called, \code{x} has to
    have two columns.}
  \item{w}{ if \code{An.biv} is called, a vector of points in [0,1]
    where to evaluate the estimated bivariate Pickands dependence
    function. If the multivariate estimator \code{An} is used instead,
    \code{w} needs to be a matrix with the same number of columns as
    \code{x} whose lines are elements of the multivariate unit simplex
    (see the last reference).}
  \item{estimator}{specifies which nonparametric rank-based estimator
    of the unknown Pickands dependence function to use in the bivariate
    case; can be either
    \code{"CFG"}(Cap\enc{é}{e}ra\enc{à}{a}-Foug\enc{è}{e}res-Genest) or
    \code{"Pickands"}.}
  \item{corrected}{ TRUE means that the bivariate estimators will be
    corrected to ensure that their value at 0 and 1 is 1.}
}
\details{
  More details can be found in the references.
}
\value{
  \code{An.biv()} returns a vector containing the values of the
  estimated Pickands dependence function at the points in \code{w} (and
  is the same as former \code{Anfun()}).

  The function \code{An} computes simultaneously the three corrected
  multivariate estimators studied in Gudendorf and Segers (2011) at the
  points in \code{w} and retuns a list whose components are
  \item{P}{ values of the Pickands estimator at the points in \code{w}.}
  \item{CFG}{ values of the CFG estimator at the points in \code{w}.}
  \item{HT}{ values of the Hall-Tajvidi estimator at the points in \code{w}.}
}
\references{
  C. Genest and J. Segers (2009).
  Rank-based inference for bivariate extreme-value copulas.
  \emph{Annals of Statistics} \bold{37}, 2990--3022.

  G. Gudendorf and J. Segers (2011). Nonparametric estimation of
  multivariate extreme-value copulas. \emph{arXiv:1107.2410v1}.
}
\seealso{\code{\link{evCopula}}, \code{\link{A}}, and \code{\link{evTestA}}.
  Further, \code{\link{evTestC}}, \code{\link{evTestK}},
  \code{\link{exchEVTest}}, and \code{\link{gofEVCopula}}. }
\examples{
## True Pickands dependence functions
curve(A(gumbelCopula(4   ), x), 0, 1)
curve(A(gumbelCopula(2   ), x), add=TRUE, col=2)
curve(A(gumbelCopula(1.33), x), add=TRUE, col=3)

## CFG estimator
curve(An.biv(rCopula(1000, gumbelCopula(4   )), x), lty=2, add=TRUE)
curve(An.biv(rCopula(1000, gumbelCopula(2   )), x), lty=2, add=TRUE, col=2)
curve(An.biv(rCopula(1000, gumbelCopula(1.33)), x), lty=2, add=TRUE, col=3)

## Pickands estimator
curve(An.biv(rCopula(1000, gumbelCopula(4   )), x, estimator="Pickands"),
      lty=3, add=TRUE)
curve(An.biv(rCopula(1000, gumbelCopula(2   )), x, estimator="Pickands"),
      lty=3, add=TRUE, col=2)
curve(An.biv(rCopula(1000, gumbelCopula(1.33)), x, estimator="Pickands"),
      lty=3, add=TRUE, col=3)

## Relationship between An.biv and An
u <- c(runif(100),0,1) # include 0 and 1
x <- rCopula(1000, gumbelCopula(4))
r <- An(x, cbind(1-u, u))
all.equal(r$CFG, An.biv(x, u))
all.equal(r$P, An.biv(x, u, estimator="Pickands"))

## A trivariate example
x <- rCopula(1000, gumbelCopula(4, dim = 3))
u <- matrix(runif(300), 100, 3)
w <- u / apply(u, 1, sum)
r <- An(x, w)

## Endpoint corrections are applied
An(x, cbind(1, 0, 0))
An(x, cbind(0, 1, 0))
An(x, cbind(0, 0, 1))

}
\keyword{htest}
\keyword{multivariate}
