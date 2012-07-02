\name{Cn}
\alias{Cn}
\title{Computes the empirical copula}
\description{
  Given a random sample from a distribution with continuous margins and
  copula C, the empirical copula is a natural nonparametric estimator of
  C.
}
\usage{
Cn(x, w)
}
\arguments{
  \item{x}{ a data matrix whose lines represent a multivariate random
    sample that will be transformed to pseudo-observations. }
  \item{w}{ a matrix with the same number of columns as \code{x} whose
    lines represent the points where the empirical copula of \code{x}
    will be evaluated.}
}
\details{
  There are several asymptotically equivalent definitions of the
  empirical copula.  Here, the empirical copula is simply defined as the
  empirical cumulative distribution function computed from the
  pseudo-observations, i.e., multivariate ranks divided by the number
  of observations plus one.
}
\value{
  Returns a \code{numeric} containing the values of the empirical copula
  of \code{x} at the points specified in \code{w}.
}
\references{
  R\enc{ü}{u}schendorf, L. (1976). Asymptotic distributions of
  multivariate rank order statistics, \emph{Annals of
    Statistics}. \bold{4}, 912--923.

  Deheuvels, P. (1979). La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, \emph{Acad. Roy.
    Belg. Bull. Cl. Sci.}, 5th Ser. \bold{65}, 274--292.

  Deheuvels, P. (1981), A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris}. \bold{26}, 29--50.
}
\seealso{\code{\link{pCopula}}.}
\examples{
n <- 20
gumbel.cop <- gumbelCopula(4)
u <- rCopula(n, gumbel.cop)
## points were to evaluate the empirical copula
w <- matrix(runif(2 * n), n, 2)
ec <- Cn(u, w)
## compare with true c.d.f.
mean(abs(pCopula(w, gumbel.cop) - ec))
## increase n to decrease the mean absolute error

## compare the empirical copula and the true copula c.d.f.
## on the diagonal of the unit square
curve(Cn(u, cbind(x, x)), 0, 1)
curve(pCopula(cbind(x, x), gumbel.cop), lty = 2, add = TRUE)
}
\keyword{multivariate}
