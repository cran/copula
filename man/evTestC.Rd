\name{evTestC}
\alias{evTestC}
\title{Large-sample Test of Multivariate Extreme-Value Dependence}
\description{
  Test of multivariate extreme-value dependence based on the empirical
  copula and max-stability. The test statistics are defined in the second
  reference. Approximate p-values for the test statistics are obtained
  by means of a \emph{multiplier} technique.
}
\usage{
evTestC(x, N = 1000)
}
\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{ number of multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
}
\details{
  More details are available in the second reference.
  See also Remillard and Scaillet (2009).
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
}
\references{
  R\enc{é}{e}millard, B. and Scaillet, O. (2009). Testing for equality
  between two copulas. \emph{Journal of Multivariate Analysis}, 100(3),
  pages 377-386.

  Kojadinovic, I., Segers, J., and Yan, J. (2011). Large-sample tests of
  extreme-value dependence for multivariate copulas. \emph{The Canadian
    Journal of Statistics} \bold{39}, 4, pages 703-720.

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
\seealso{\code{\link{evTestK}}, \code{\link{evTestA}}, \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}.}
\examples{
## Do these data come from an extreme-value copula?
evTestC(rCopula(200, gumbelCopula(3)))
evTestC(rCopula(200, claytonCopula(3)))

## Three-dimensional examples
evTestC(rCopula(200, gumbelCopula(3, dim=3)))
evTestC(rCopula(200, claytonCopula(3, dim=3)))
\dontshow{
 set.seed(101)
 G.t <- evTestC(rCopula(200, gumbelCopula(3, dim=3)))
 C.t <- evTestC(rCopula(200, claytonCopula(3, dim=3)))
 eT3 <- evTestC(rCopula(200, tevCopula(.8, df=3)))
 stopifnot(all.equal(G.t$p.value, 0.1543, tol=.001),
           all.equal(C.t$p.value, 4995/9999990, tol= 1e-7),
           all.equal(eT3$p.value, 0.407092907092907, tol= 1e-7))
}
}
\keyword{htest}
\keyword{multivariate}
