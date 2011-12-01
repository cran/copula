\name{exchEVTest}
\alias{exchEVTest}

\title{Test of exchangeability for extreme-value and left-tail
  decreasing bivariate copulas}

\description{
  Test for assessing the exchangeability of the underlying
  bivariate copula when it is either extreme-value or left-tail
  decreasing. The test uses the nonparametric estimators of the Pickands
  dependence function studied by Genest and Segers (2009).
  The test statistic is defined in the second reference.
  An approximate p-value for the test statistic is obtained
  by means of a \emph{multiplier} technique.  
}
\usage{
exchEVTest(x, N = 1000,  estimator = "CFG", derivatives = "Cn", m = 100)
}

\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations.}
  \item{N}{ number of multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{estimator}{ specifies which nonparametric estimator of the Pickands
    dependence function to use; can be either \code{"CFG"} or
    \code{"Pickands"}; see Genest and Segers (2009).}
  \item{derivatives}{ specifies how the derivatives of the unknown 
    copula are estimated; can be either \code{"An"} or \code{"Cn"}.
    The former should be used under the assumption of extreme-value
    dependence. The latter is faster. See the second reference.}
  \item{m}{ specifies the size of the integration grid for the statistic.}
}

\details{

  More details are available in the two first references.
  
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{pvalue}{ corresponding approximate p-value. }
}
\references{

  C. Genest and J. Segers (2009). Rank-based inference for bivariate
  extreme-value copulas. \emph{Annals of Statistics}, 37, pages 2990-3022.
  
  I. Kojadinovic and J. Yan (2012). A nonparametric test of
  exchangeability for extreme-value and left-tail
  decreasing bivariate copulas. \emph{The Scandinavian Journal of
    Statistics}. In press.

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
  
\seealso{ \code{\link{exchTest}}, \code{\link{gofCopula}}. }
\examples{
## Do these data come from exchangeable copulas? 
exchEVTest(rcopula(gumbelCopula(3), 200)) 
exchEVTest(rcopula(claytonCopula(3), 200))

## Creating asymmetric data 
khoudraji <- function(cop,n,a=0.6,b=0.95)
  {
    u <- rcopula(cop,n)
    v <- matrix(runif(2*n),n,2)
    x <- cbind(pmax(u[,1]^(1/a),v[,1]^(1/(1-a))),
               pmax(u[,2]^(1/b),v[,2]^(1/(1-b))))
    x
  }
exchEVTest(khoudraji(gumbelCopula(3),200))
exchEVTest(khoudraji(claytonCopula(3),200))
}

\keyword{htest}
\keyword{multivariate}
