\name{radSymTest}
\alias{radSymTest}
\title{Test of Exchangeability for a Bivariate Copula}
\description{

  Test for assessing the radial symmetry of the underlying multivariate
  copula based on the empirical copula.  The test statistic is a
  multivariate extension of the definition adopted in the first
  reference. An approximate p-value for the test statistic is obtained
  by means of a appropriate \emph{bootstrap} which can take the presence
  of ties in the component series of the data into accont; see the
  second reference.

}
\usage{
radSymTest(x, N = 1000, ties = NA)
}
\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations.}
  \item{N}{ number of boostrap iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{ties}{ logical; if \code{TRUE}, the boostrap procedure is
    adapted to the presence of ties in any of the coordinate samples
    of \code{x}; the default value of \code{NA} indicates that the
    presence/absence of ties will be checked for automatically.}
}
\details{
  More details are available in the second reference.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{

  Genest, C. and G. \enc{Nešlehová}{Neslehova}, J. (2014).  On tests of
  radial symmetry for bivariate copulas. \emph{Statistical Papers}
  \bold{55}, 1107--1119.

  Kojadinovic, I. (2016). Some copula inference procedures adapted to
  the presence of ties. \url{http://arxiv.org/abs/1609.05519}.

}
\seealso{ \code{\link{exchTest}}, \code{\link{exchEVTest}}, \code{\link{gofCopula}}. }
\examples{
## Data from radially symmetric copulas
radSymTest(rCopula(200, frankCopula(3)))
\donttest{radSymTest(rCopula(200, normalCopula(0.7, dim = 3)))}

## Data from non radially symmetric copulas
radSymTest(rCopula(200, claytonCopula(3)))
\donttest{radSymTest(rCopula(200, gumbelCopula(2, dim=3)))}
}
\keyword{htest}
\keyword{multivariate}
