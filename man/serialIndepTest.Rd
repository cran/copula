\name{serialIndepTest}
\alias{serialIndepTest}
\alias{serialIndepTestSim}
\title{Serial Independence Test for Continuous Time Series Based
  on the Empirical Copula Process}
\description{
  Serial independence test based on the empirical copula process as
  proposed in Ghoudi et al. (2001) and Genest and Rémillard (2004). The
  test, which is the serial analog of \code{\link{indepTest}}, can be
  seen as composed of three steps: (i) a simulation step, which consists
  in simulating the distribution of the test statistics under serial
  independence for the sample size under consideration; (ii) the test
  itself, which consists in computing the approximate p-values of the
  test statistics with respect to the empirical distributions obtained
  in step (i); and (iii) the display of a graphic, called a
  \emph{dependogram}, enabling to understand the type of departure from
  serial independence, if any. More details can be found in the articles
  cited in the reference section.
}
\usage{
  serialIndepTestSim(n, lag.max, m=lag.max+1, N=1000, verbose = TRUE,
                     print.every = NULL)
  serialIndepTest(x, d, alpha=0.05)
}
\arguments{
  \item{n}{ length of the time series when simulating the distribution of
    the test statistics under serial independence.}
  \item{lag.max}{ maximum lag.}
  \item{m}{ maximum cardinality of the subsets of 'lags' for which a test
    statistic is to be computed. It makes sense to consider \code{m <<
      lag.max+1} especially when \code{lag.max} is large.}
  \item{N}{ number of repetitions when simulating under serial
    independence.}
   \item{print.every}{ is deprecated in favor of \code{verbose}.}
   \item{verbose}{ a logical specifying if progress
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{x}{ numeric vector containing the time series whose serial
    independence is to be tested.}
  \item{d}{ object of class \code{serialIndepTestDist} as returned by the
    function \code{serialIndepTestSim}. It can be regarded as the
    empirical distribution of the test statistics under serial
    independence.}
  \item{alpha}{ significance level used in the computation of the
    critical values for the test statistics.}
}
\details{See the references below for more details, especially the third
  and fourth ones.}
\value{The function \code{"serialIndepTestSim"} returns an object of class
   \code{"serialIndepTestDist"} whose attributes are: \code{sample.size},
  \code{lag.max}, \code{max.card.subsets}, \code{number.repetitons},
  \code{subsets} (list of the subsets for which test statistics have
  been computed), \code{subsets.binary} (subsets in binary 'integer'
  notation), \code{dist.statistics.independence} (a \code{N} line matrix
  containing the values of the test statistics for each subset and each
  repetition) and \code{dist.global.statistic.independence} (a vector a
  length \code{N} containing the values of the serial version of the
  global Cramér-von Mises test statistic for each repetition --- see
  last reference p.175).

  The function \code{"serialIndepTest"} returns an object of class
  \code{"indepTest"} whose attributes are: \code{subsets},
  \code{statistics}, \code{critical.values}, \code{pvalues},
  \code{fisher.pvalue} (a p-value resulting from a combination \emph{à la}
  Fisher of the subset statistic p-values), \code{tippett.pvalue} (a
  p-value resulting from a combination \emph{à la} Tippett of the subset
  statistic p-values), \code{alpha} (global significance level of the
  test), \code{beta} (1 - \code{beta} is the significance level per
  statistic), \code{global.statistic} (value of the global Cramér-von
  Mises statistic derived directly from the serial independence
  empirical copula process --- see last reference p 175) and
  \code{global.statistic.pvalue} (corresponding p-value).
}
\references{
  Deheuvels, P. (1979). La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, \emph{Acad. Roy.
    Belg. Bull. Cl. Sci.}, 5th Ser. 65:274--292.

  Deheuvels, P. (1981), A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris}. 26:29--50.

  Genest, C. and R\enc{é}{e}millard, B. (2004), Tests of independence and
  randomness based on the empirical copula process. \emph{Test}
  \bold{13}, 335--369.

  Genest, C., Quessy, J.-F., and R\enc{é}{e}millard, B. (2006). Local
  efficiency of a Cramer-von Mises test of independence,
  \emph{Journal of Multivariate Analysis} \bold{97}, 274--294.

  Genest, C., Quessy, J.-F., and R\enc{é}{e}millard, B. (2007), Asymptotic
  local efficiency of Cramér-von Mises tests for multivariate
  independence. \emph{The Annals of Statistics} \bold{35}, 166--191.
}
\seealso{
  \code{\link{indepTest}},
  \code{\link{multIndepTest}},
  \code{\link{multSerialIndepTest}},
  \code{\link{dependogram}}
}
\examples{
## AR 1 process

ar <-  numeric(200)
ar[1] <- rnorm(1)
for (i in 2:200)
  ar[i] <- 0.5 * ar[i-1] + rnorm(1)
x <- ar[101:200]

## In order to test for serial independence, the first step consists
## in simulating the distribution of the test statistics under
## serial independence for the same sample size, i.e. n=100.
## As we are going to consider lags up to 3, i.e., subsets of
## {1,...,4} whose cardinality is between 2 and 4 containing {1},
## we set lag.max=3. This may take a while...

d <- serialIndepTestSim(100,3)

## The next step consists in performing the test itself:

test <- serialIndepTest(x,d)

## Let us see the results:

test

## Display the dependogram:

dependogram(test,print=TRUE)

## NB: In order to save d for future use, the save function can be used.
}
\keyword{htest}
