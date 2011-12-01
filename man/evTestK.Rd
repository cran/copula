\name{evTestK}
\alias{evTestK}

\title{Bivariate test of extreme-value dependence based on Kendall's process}
\description{
  Test of extreme-value dependence based on the  bivariate probability integral
  transformation. The test statistic is defined in
  Ben Ghorbal, Nešlehová and Genest (2009). 
}
\usage{
evTestK(x, method = "fsample")
}

\arguments{
  \item{x}{ a data matrix. }
  \item{method}{ specifies the variance estimation method; can be either
    \code{"fsample"} (finite-sample), \code{"asymptotic"} or \code{"jackknife"}.}
}
\details{
  The code for this test was generously provided by Johanna Nešlehová. 
  More details are available in Appendix B of
  Ben Ghorbal, Nešlehová and Genest (2009).
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{pvalue}{ corresponding p-value. }
}
\references{
  M. Ben Ghorbal, C. Genest, and J. Nešlehová (2009). On the test of Ghoudi,
  Khoudraji, and Rivest for extreme-value dependence. \emph{The Canadian
    Journal of Statistics}, 37, pages 1-9.

}
\seealso{ \code{\link{evTestC}}, \code{\link{evTestA}}, \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{Anfun}}. }
\examples{
## Do the data come from an extreme-value copula? 
evTestK(rcopula(gumbelCopula(3), 200)) 
evTestK(rcopula(claytonCopula(3), 200))
}

\encoding{UTF-8}
\keyword{htest}
\keyword{multivariate}
