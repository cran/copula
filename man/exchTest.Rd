\name{exchTest}
\alias{exchTest}
\title{Test of Exchangeability for a Bivariate Copula}
\description{
  Test for assessing the exchangeability of the underlying
  bivariate copula based on the empirical copula.
  The test statistics are defined in the first two references.
  Approximate p-values for the test statistics are obtained
  by means of a \emph{multiplier} technique if there are no ties in the
  component series of the bivariate data, or by means of an appropriate
  bootstrap otherwise.
}
\usage{
exchTest(x, N = 1000, ties = NA,
         ties.method = eval(formals(rank)$ties.method), m = 0)
}
\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations.}
  \item{N}{ number of multiplier or boostrap iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{ties}{ logical; if \code{FALSE}, approximate p-values are computed by
    means of a multiplier bootstrap; if \code{TRUE}, a boostrap adapted to the
    presence of ties in any of the coordinate samples of \code{x}
    is used; the default value of \code{NA} indicates that the
    presence/absence of ties will be checked for automatically.}
  \item{ties.method}{string specifying how ranks should be computed if
    there are ties in any of the coordinate samples of \code{x}; passed
    to \code{\link{pobs}}.}
  \item{m}{ if \code{m=0}, integration in the Cramér--von Mises
    statistic is carried out with respect to the empirical copula; if
  \code{m > 0}, integration is carried out with respect to the Lebesgue
  measure and \code{m} specifies the size of the integration grid.}
}
\details{
  More details are available in the references.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Genest, C., G. \enc{Nešlehová}{Neslehova}, J. and Quessy,
  J.-F. (2012).  Tests of symmetry for bivariate copulas. \emph{Annals
  of the Institute of Statistical Mathematics} \bold{64}, 811--834.

  Kojadinovic, I. and Yan, J. (2012). A nonparametric test of
  exchangeability for extreme-value and left-tail
  decreasing bivariate copulas. \emph{The Scandinavian Journal of
    Statistics} \bold{39:3}, 480--496.

  Kojadinovic, I. (2017). Some copula inference procedures adapted to
  the presence of ties. \emph{Computational Statistics and Data
    Analysis} \bold{112}, 24--41, \url{http://arxiv.org/abs/1609.05519}.
}
\seealso{ \code{\link{radSymTest}}, \code{\link{exchEVTest}}, \code{\link{gofCopula}}. }
\examples{
## Data from an exchangeable copulas
exchTest(rCopula(200,  gumbelCopula(3)))
\donttest{exchTest(rCopula(200, claytonCopula(3)))}

## An asymmetric Khoudraji-Clayton copula
kc <- khoudrajiCopula(copula1 = indepCopula(),
                      copula2 = claytonCopula(6),
                      shapes = c(0.4, 0.95))
exchTest(rCopula(200, kc))
}
\keyword{htest}
\keyword{multivariate}
