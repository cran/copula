\name{indepTest}
\alias{indepTest}
\alias{indepTestSim}
\alias{dependogram}

\title{Independence test among continuous random variables based
  on the empirical copula process}

\description{
  Multivariate independence test based on the empirical 
  copula process as proposed by Christian Genest and Bruno 
  Rémillard. The test can be seen as composed of
  three steps: (i) a simulation step, which consists of simulating the
  distribution of the test statistics under independence for the sample
  size under consideration; (ii) the test itself, which consists of
  computing the approximate p-values of the test statistics with respect
  to the empirical distributions obtained in step (i); and (iii) the
  display of a graphic, called a \emph{dependogram}, enabling to
  understand the type of departure from independence, if any. More details can
  be found in the articles cited in the reference section.
}

\usage{
indepTestSim(n, p, m = p, N = 1000, print.every = 100)
indepTest(x, d, alpha=0.05)
dependogram(test, pvalues = FALSE, print = FALSE)
}

\arguments{
  \item{n}{Sample size when simulating the distribution of the test
    statistics under independence.}
  \item{p}{Dimension of the data when simulating the distribution of the test
    statistics under independence.}
  \item{m}{Maximum cardinality of the subsets of variables for which a test statistic
    is to be computed. It makes sense to consider \code{m << p} especially when
    \code{p} is large.}
  \item{N}{Number of repetitions when simulating under independence.}
  \item{print.every}{Progress is printed every \code{"print.every"}
    iterations. No progress is printed if it is nonpositive.}
  \item{x}{Data frame or data matrix containing realizations (one per
    line) of the random vector whose independence is to be tested.}
  \item{d}{Object of class \code{"indepTestDist"} as returned by
  the function \code{"indepTestSim"}. It can be regarded as the empirical distribution
  of the test statistics under independence.}
  \item{alpha}{Significance level used in the computation of the critical values for
    the test statistics.}
  \item{test}{Object of class \code{"indepTest"} as return by the function \code{indepTest}.}
  \item{pvalues}{Logical indicating whether the dependogram should be
    drew from test statistics or the corresponding p-values.}
  \item{print}{Logical indicating whether details should be printed.}
}
\details{See the references below for more details, especially the third
  one.
}

\value{The function \code{"indepTestSim"} returns an object of class
  \code{"indepTestDist"} whose attributes are: \code{sample.size},
  \code{data.dimension}, \code{max.card.subsets},
  \code{number.repetitons}, \code{subsets} (list of the subsets for
  which test statistics have been computed), \code{subsets.binary}
  (subsets in binary 'integer' notation),
  \code{dist.statistics.independence} (a \code{N} line matrix containing
  the values of the test statistics for each subset and each repetition)
  and \code{dist.global.statistic.independence} (a vector a length
  \code{N} containing the values of the global Cramér-von Mises test
  statistic for each repetition - see last reference p. 175).
  
  The function \code{"indepTest"} returns an object of class
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
  process - see last reference p. 175) and
  \code{global.statistic.pvalue} (corresponding p-value).
}

\references{
  P. Deheuvels (1979). La fonction de dépendance empirique et ses
  propriétés: un test non paramétrique d'indépendance, \emph{Acad. Roy.
    Belg. Bull. Cl. Sci.}, 5th Ser. 65:274--292.
  
  P. Deheuvels (1981), A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris}. 26:29--50.
  
  C. Genest and B. Rémillard (2004), Tests of independence and
  randomness based on the empirical copula process. \emph{Test}, 13:335--369.
  
  C. Genest, J.-F. Quessy and B. Rémillard (2006). Local
  efficiency of a Cramer-von Mises test of independence,
  \emph{Journal of Multivariate Analysis}, 97:274--294.
  
  C. Genest, J.-F. Quessy and B. Rémillard (2007), Asymptotic
  local efficiency of Cramér-von Mises tests for multivariate
  independence. \emph{The Annals of Statistics}, 35:166--191.
}

\seealso{
  \code{\link{serialIndepTest}}, 
  \code{\link{multIndepTest}}, 
  \code{\link{multSerialIndepTest}}.
}
\examples{
## Consider the following example taken from
## Genest and Remillard (2004), p 352:

x <- matrix(rnorm(500),100,5)
x[,1] <- abs(x[,1]) * sign(x[,2] * x[,3])
x[,5] <- x[,4]/2 + sqrt(3) * x[,5]/2

## In order to test for independence "within" x, the first step consists
## in simulating the distribution of the test statistics under
## independence for the same sample size and dimension,
## i.e. n=100 and p=5. As we are going to consider all the subsets of
## {1,...,5} whose cardinality is between 2 and 5, we set p=m=5.
## This may take a while...

d <- indepTestSim(100,5)

## The next step consists of performing the test itself:

test <- indepTest(x,d)

## Let us see the results:

test

## Display the dependogram with the details:

dependogram(test,print=TRUE)

## We could have tested for a weaker form of independence, for instance,
## by only computing statistics for subsets whose cardinality is between 2
## and 3. Consider for instance the following data:
y <- matrix(runif(500),100,5)
## and perform the test:
d <- indepTestSim(100,5,3)
test <- indepTest(y,d)
test
dependogram(test,print=TRUE)

## NB: In order to save d for future use, the save function can be used.
}

\encoding{UTF-8}
\keyword{htest}
