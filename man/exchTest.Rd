\name{exchTest}
\alias{exchTest}
\title{Test of Exchangeability for a Bivariate Copula}
\description{
  Test for assessing the exchangeability of the underlying
  bivariate copula based on the empirical copula.
  The test statistics are defined in the first two references.
  Approximate p-values for the test statistics are obtained
  by means of a \emph{multiplier} technique.
}
\usage{
exchTest(x, N = 1000, m = 0)
}
\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations.}
  \item{N}{ number of multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{m}{ if \code{m=0}, integration in the Cramér--von Mises
    statistic is carried out with respect to the empirical copula; if
  \code{m > 0}, integration is carried out with respect to the Lebesgue
  measure and \code{m} specifies the size of the integration grid.}
}
\details{
  More details are available in the two first references.
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
}
\references{
  Genest, C., \enc{Nešlehová}{Neslehova}, J. and Quessy, J.-F. (2012). Tests of symmetry for
  bivariate copulas. \emph{Annals of the Institute of Statistical
  Mathematics}. In press.

  Kojadinovic, I. and Yan, J. (2012). A nonparametric test of
  exchangeability for extreme-value and left-tail
  decreasing bivariate copulas. \emph{The Scandinavian Journal of
    Statistics}. In press.

  Kojadinovic, I. and Yan, J. (2010). Modeling Multivariate Distributions
  with Continuous Margins Using the copula R Package. \emph{Journal of
  Statistical Software}, 34(9), pages 1-20.
}
\note{
  This test was derived under the assumption of continuous margins,
  which implies that ties occur with probability zero. The
  presence of ties in the data might substantially affect the
  approximate p-value. One way of dealing with ties was suggested in the
  last reference.
}
\seealso{ \code{\link{exchEVTest}}, \code{\link{gofCopula}}. }
\examples{
## Do these data come from exchangeable copulas?
exchTest(rCopula(200, gumbelCopula(3)))
exchTest(rCopula(200, claytonCopula(3)))

## Creating asymmetric data
khoudraji <- function(cop,n,a=0.6,b=0.95)
  {
    u <- rCopula(n, cop)
    v <- matrix(runif(2*n),n,2)
    x <- cbind(pmax(u[,1]^(1/a),v[,1]^(1/(1-a))),
               pmax(u[,2]^(1/b),v[,2]^(1/(1-b))))
    x
  }
exchTest(khoudraji(gumbelCopula(3),200))
exchTest(khoudraji(claytonCopula(3),200))
}
\keyword{htest}
\keyword{multivariate}
