\name{multSerialIndepTest}
\alias{multSerialIndepTest}
\title{Serial Independence Test for Multivariate Time Series via Empirical Copula}
\description{
  Analog of the serial independence test based on the empirical
  copula process proposed by Christian Genest and Bruno Rémillard (see
  \code{\link{serialIndepTest}}) for \emph{multivariate} time
  series.  The main difference comes from the fact that critical values
  and p-values are obtained through the bootstrap/permutation
  methodology, since, here, test statistics are not distribution-free.
}
\usage{
multSerialIndepTest(x, lag.max, m = lag.max+1, N = 1000, alpha = 0.05,
                    verbose = interactive(), print.every = NULL)
}
\arguments{
  \item{x}{data frame or matrix of multivariate continuous time series
    whose serial independence is to be tested.}
  \item{lag.max}{ maximum lag.}
  \item{m}{ maximum cardinality of the subsets of 'lags' for
    which a test statistic is to be computed.  It makes sense to consider
    \code{m << lag.max+1} especially when \code{lag.max} is large.}
  \item{N}{ number of bootstrap/permutation samples.}
  \item{alpha}{ significance level used in the computation of the
    critical values for the test statistics.}
  \item{print.every}{ is deprecated in favor of \code{verbose}.}
  \item{verbose}{ a logical specifying if progress
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
}
\details{See the references below for more details, especially the last
  one.
}
\value{
  The function \code{"multSerialIndepTest"} returns an object of class
  \code{"indepTest"} whose attributes are: \code{subsets},
  \code{statistics}, \code{critical.values}, \code{pvalues},
  \code{fisher.pvalue} (a p-value resulting from a combination \emph{à la}
    Fisher of the subset statistic p-values), \code{tippett.pvalue} (a p-value
  resulting from a combination \emph{à la} Tippett of the subset
  statistic p-values), \code{alpha} (global significance level of the
  test), \code{beta} (1 - \code{beta} is the significance level per statistic),
  \code{global.statistic} (value of the global Cramér-von Mises
  statistic derived directly from
  the independence empirical copula process - see \code{In} in the last
  reference) and \code{global.statistic.pvalue} (corresponding p-value).
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Deheuvels, P. (1979)  La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance.
  \emph{Acad. Roy. Belg. Bull. Cl. Sci.}, 5th Ser. \bold{65}, 274--292.

  Deheuvels, P. (1981)  A non parametric test for independence.
  \emph{Publ. Inst. Statist. Univ. Paris} \bold{26}, 29--50.

  Genest, C. and R\enc{é}{e}millard, B. (2004)
  Tests of independence and randomness based on the empirical copula
  process.
  \emph{Test} \bold{13}, 335--369.

  Ghoudi, K., Kulperger, R., and R\enc{é}{e}millard, B. (2001)
  A nonparametric test of serial independence for times series and residuals.
  \emph{Journal of Multivariate Analysis} \bold{79}, 191--218.

  Kojadinovic, I. and Yan, J. (2011)
  Tests of multivariate serial independence based on a Möbius
  decomposition of the independence empirical copula process.
  \emph{Annals of the Institute of Statistical Mathematics} \bold{63}, 347--373.
}
\seealso{
  \code{\link{serialIndepTest}},
  \code{\link{indepTest}},
  \code{\link{multIndepTest}},
  \code{\link{dependogram}}
}
\examples{
## A multivariate time series {minimal example for demo purposes}
d <- 2
n <- 100 # sample size *and* "burn-in" size
param <- 0.25
A <- matrix(param,d,d) # the bivariate AR(1)-matrix
set.seed(17)
ar <- matrix(rnorm(2*n * d), 2*n,d) # used as innovations
for (i in 2:(2*n))
  ar[i,] <- A \%*\% ar[i-1,] + ar[i,]
## drop burn-in :
x <- ar[(n+1):(2*n),]

## Run the test
test <- multSerialIndepTest(x,3)
test

## Display the dependogram
dependogram(test,print=TRUE)
}
\keyword{htest}
