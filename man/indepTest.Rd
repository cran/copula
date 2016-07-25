\name{indepTest}
\alias{indepTest}
\alias{indepTestSim}
\alias{dependogram}
\title{Test Independence of Continuous Random Variables via Empirical Copula}
\description{
  Multivariate independence test based on the empirical
  copula process as proposed by Christian Genest and Bruno
  R\enc{é}{e}millard.  The test can be seen as composed of
  three steps: (i) a simulation step, which consists of simulating the
  distribution of the test statistics under independence for the sample
  size under consideration; (ii) the test itself, which consists of
  computing the approximate p-values of the test statistics with respect
  to the empirical distributions obtained in step (i); and (iii) the
  display of a graphic, called a \emph{dependogram}, enabling to
  understand the type of departure from independence, if any.  More
  details can be found in the articles cited in the reference section.
}
\usage{
indepTestSim(n, p, m = p, N = 1000, verbose = interactive(), print.every = NULL)
indepTest(x, d, alpha=0.05)
dependogram(test, pvalues = FALSE, print = FALSE)
}
\arguments{
%% indepTestSim():
  \item{n}{ sample size when simulating the distribution of the test
    statistics under independence.}
  \item{p}{ dimension of the data when simulating the distribution of the test
    statistics under independence.}
  \item{m}{maximum cardinality of the subsets of variables for which a
    test statistic is to be computed.  It makes sense to consider \eqn{m
      \ll p}{m << p} especially when \code{p} is large.}
  \item{N}{number of repetitions when simulating under independence.}
  \item{print.every}{ is deprecated in favor of \code{verbose}.}
  \item{verbose}{a logical specifying if progress
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
%% indepTest():
  \item{x}{ data frame or data matrix containing realizations (one per
    line) of the random vector whose independence is to be tested.}
  \item{d}{ object of class \code{"indepTestDist"} as returned by
    the function \code{indepTestSim()}.  It can be regarded as the
    empirical distribution of the test statistics under independence.}
  \item{alpha}{ significance level used in the computation of the
    critical values for the test statistics.}
%% dependogram():
  \item{test}{object of class \code{"indepTest"} as returned by \code{indepTest()}.}
  \item{pvalues}{ logical indicating whether the dependogram should be
    drew from test statistics or the corresponding p-values.}
  \item{print}{ logical indicating whether details should be printed.}
}
\details{
  The current (C code) implementation of \code{indepTestSim()}
  uses (RAM) memory of size \eqn{O(n^2 p)}{O(n^2 * p)}, and time
  % ../src/ecIndepTest.c : J_size, in simulate_empirical_copula()
  \eqn{O(N n^2 p)}{O(N * n^2 * p)}.  This renders it unfeasible when
  \code{n} is large.% "vague" on purpose; n=10'000 is barely "ok" w/ 8GB RAM and > 1 hour waiting

  See the references below for more details, especially Genest
  and R\enc{é}{e}millard (2004).
}
\value{The function \code{indepTestSim()} returns an object of class
  \code{"indepTestDist"} whose attributes are: \code{sample.size},
  \code{data.dimension}, \code{max.card.subsets},
  \code{number.repetitons}, \code{subsets} (list of the subsets for
  which test statistics have been computed), \code{subsets.binary}
  (subsets in binary 'integer' notation),
  \code{dist.statistics.independence} (a \code{N} line matrix containing
  the values of the test statistics for each subset and each repetition)
  and \code{dist.global.statistic.independence} (a vector a length
  \code{N} containing the values of the global Cramér-von Mises test
  statistic for each repetition -- see Genest \emph{et al} (2007), p.175).

  The function \code{indepTest()} returns an object of class
  \code{"indepTest"} whose attributes are: \code{subsets},
  \code{statistics}, \code{critical.values}, \code{pvalues},
  \code{fisher.pvalue} (a p-value resulting from a combination \emph{à la}
  Fisher of the subset statistic p-values), \code{tippett.pvalue} (a p-value
  resulting from a combination \emph{à la} Tippett of the subset
  statistic p-values),
  \code{alpha} (global significance level of the test), \code{beta}
  (1 - \code{beta} is the significance level per statistic),
  \code{global.statistic} (value of the global Cramér-von Mises
  statistic derived directly from the independence empirical copula
  process - see Genest \emph{et al} (2007), p.175) and
  \code{global.statistic.pvalue} (corresponding p-value).
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{%% do not add or remove easily
  Deheuvels, P. (1979).
  La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance,
  \emph{Acad. Roy. Belg. Bull. Cl. Sci.}, 5th Ser. \bold{65}, 274--292.

  Deheuvels, P. (1981)
  A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris}. \bold{26}, 29--50.

  Genest, C. and R\enc{é}{e}millard, B. (2004)
  Tests of independence and randomness based on the empirical copula process.
  \emph{Test} \bold{13}, 335--369.

  Genest, C., Quessy, J.-F., and R\enc{é}{e}millard, B. (2006).
  Local efficiency of a Cramer-von Mises test of independence,
  \emph{Journal of Multivariate Analysis} \bold{97}, 274--294.

  Genest, C., Quessy, J.-F., and R\enc{é}{e}millard, B. (2007)
  Asymptotic local efficiency of Cramér-von Mises tests for multivariate independence.
  \emph{The Annals of Statistics} \bold{35}, 166--191.
}
\seealso{
  \code{\link{serialIndepTest}},
  \code{\link{multIndepTest}},
  \code{\link{multSerialIndepTest}}.
}
\examples{
## Consider the following example taken from
## Genest and Remillard (2004), p 352:

set.seed(2004)
x <- matrix(rnorm(500),100,5)
x[,1] <- abs(x[,1]) * sign(x[,2] * x[,3])
x[,5] <- x[,4]/2 + sqrt(3) * x[,5]/2

## In order to test for independence "within" x, the first step consists
## in simulating the distribution of the test statistics under
## independence for the same sample size and dimension,
## i.e. n=100 and p=5. As we are going to consider all the subsets of
## {1,...,5} whose cardinality is between 2 and 5, we set p=m=5.

## For a realistic N = 1000 (default), this takes a few seconds:
N. <- if(copula:::doExtras()) 1000 else 120
N.
system.time(d <- indepTestSim(100, 5, N = N.))
## For N=1000,  2 seconds (lynne 2015)
## You could save 'd' for future use, via  saveRDS()

## The next step consists of performing the test itself (and print its results):
(iTst <- indepTest(x,d))

## Display the dependogram with the details:
dependogram(iTst, print=TRUE)

## We could have tested for a weaker form of independence, for instance,
## by only computing statistics for subsets whose cardinality is between 2
## and 3. Consider for instance the following data:
y <- matrix(runif(500),100,5)
## and perform the test:
system.time( d <- indepTestSim(100,5,3, N=N.) )
iTy <- indepTest(y,d)
iTy
dependogram(iTy, print=TRUE)

\dontshow{
try( indepTestSim(2e9, 5) )#-> error ".. n ... too large" % segfault previously
}
}
\keyword{htest}
