\name{evTestK}
\alias{evTestK}
\title{Bivariate test of extreme-value dependence based on Kendall's process}
\description{
  Test of extreme-value dependence based on the  bivariate probability integral
  transformation. The test statistic is defined in Ben Ghorbal,
  Nešlehová and Genest (2009).
}
\usage{
evTestK(x, method = c("fsample","asymptotic","jackknife"))
}
\arguments{
  \item{x}{ a data matrix. }
  \item{method}{ specifies the variance estimation method; can be either
    \code{"fsample"} (finite-sample, the default), \code{"asymptotic"}
    or \code{"jackknife"}.}
}
\details{
  The code for this test was generously provided by Johanna Nešlehová.
  More details are available in Appendix B of Ben Ghorbal, Nešlehová and
  Genest (2009).
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding p-value. }
}
\references{
  Ghorbal, M. B., Genest, C., and Nešlehová, J. (2009)
  On the test of Ghoudi, Khoudraji, and Rivest for extreme-value dependence.
  \emph{The Canadian Journal of Statistics} \bold{37}, 1--9.
}
\seealso{ \code{\link{evTestC}}, \code{\link{evTestA}}, \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}.
}
\examples{
## Do the data come from an extreme-value copula?
evTestK(rCopula(200, gumbelCopula(3)))
evTestK(rCopula(200, claytonCopula(3)))
}
\keyword{htest}
\keyword{multivariate}
