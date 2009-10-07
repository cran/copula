\name{evTest}
\alias{evTest}

\title{Large-sample test of extreme-value dependence}
\description{
  Test of extreme-value dependence based on the empirical copula and
  max-stability. The test statistics are defined in Kojadinovic and Yan
  (2009). Approximate p-values for the test statistics are obtained
  by means of a fast \emph{multiplier} technique.  
}
\usage{
evTest(x, N = 1000)
}

\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{ number of multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
}
\details{

  More details are available in Kojadinovic and Yan (2009).
  See also Remillard and Scaillet (2009).
  
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{pvalue}{ corresponding approximate p-value. }
}
\references{

  B. Remillard and O. Scaillet (2009). Testing for equality between two
  copulas. \emph{Journal of Multivariate Analysis}, 100(3),
  pages 377-386.
  
  I. Kojadinovic and J. Yan (2009). Large-sample tests
    of extreme-value dependence for multivariate copulas.
    Submitted.
}
\seealso{ \code{\link{evCopula}}. }
\examples{
## Do the data come from an extreme-value copula? 
evTest(rcopula(gumbelCopula(3), 200)) 
evTest(rcopula(claytonCopula(3), 200))

## A three-dimensional example
evTest(rcopula(gumbelCopula(3, dim=3), 200)) 
evTest(rcopula(claytonCopula(3, dim=3), 200))
}

\keyword{htest}
\keyword{multivariate}
