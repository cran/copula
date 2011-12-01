\name{evTestC}
\alias{evTestC}

\title{Large-sample test of multivariate extreme-value dependence}
\description{
  Test of multivariate extreme-value dependence based on the empirical copula and
  max-stability. The test statistics are defined in the second reference. Approximate p-values for the test statistics are obtained
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
  \item{pvalue}{ corresponding approximate p-value. }
}
\references{

  B. Remillard and O. Scaillet (2009). Testing for equality between two
  copulas. \emph{Journal of Multivariate Analysis}, 100(3),
  pages 377-386.
  
  I. Kojadinovic, J. Segers and J. Yan (2011). Large-sample tests of
  extreme-value dependence for multivariate copulas. \emph{The Canadian
    Journal of Statistics}, 39:4, pages 703-720.
    
  I. Kojadinovic and J. Yan (2010). Modeling Multivariate Distributions
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
  
\seealso{ \code{\link{evTestK}}, \code{\link{evTestA}}, \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{Anfun}}. }
\examples{
## Do these data come from an extreme-value copula? 
evTestC(rcopula(gumbelCopula(3), 200)) 
evTestC(rcopula(claytonCopula(3), 200))

## Three-dimensional examples
evTestC(rcopula(gumbelCopula(3, dim=3), 200)) 
evTestC(rcopula(claytonCopula(3, dim=3), 200))
}

\keyword{htest}
\keyword{multivariate}
