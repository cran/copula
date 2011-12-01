\name{multSerialIndepTest}
\alias{multSerialIndepTest}


\title{Serial independence test for multivariate continuous time series
  based on the empirical copula process}

\description{
  Analog of the serial independence test based on the empirical 
  copula process proposed by Christian Genest and Bruno Rémillard (see
  \code{\link{serialIndepTest}}) for \emph{multivariate} time
  series. The main difference comes from the fact that critical values
  and p-values are obtainted through the bootstrap/permutation
  methodology, since, here, test statistics are not distribution-free.
}

\usage{
multSerialIndepTest(x, lag.max, m=lag.max+1, N=1000, alpha=0.05,
                    print.every = 100)
}

\arguments{
  \item{x}{Data frame or data matrix containing realizations the
    multivaraite continuous time series whose serial independence is to be tested.}
  \item{lag.max}{Maximum lag.}
  \item{m}{Maximum cardinality of the subsets of 'lags' for
    which a test statistic is to be computed. It makes sense to consider
    \code{m << lag.max+1} especially when \code{lag.max} is large.}
  \item{N}{Number of bootstrap/permutation samples.}
  \item{alpha}{Significance level used in the computation of the critical values for
    the test statistics.}
  \item{print.every}{Progress is printed every \code{"print.every"}
    iterations. No progress is printed if it is nonpositive.}
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
  resulting from a combination \emph{à la} Tippett of the subset statistic p-values), \code{alpha} (global significance level of the test), \code{beta}
  (1 - \code{beta} is the significance level per statistic),
  \code{global.statistic} (value of the global Cramér-von Mises
  statistic derived directly from
  the independence empirical copula process - see \code{In} in the last reference) and
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
  
  K. Ghoudi, R. Kulperger, and B. Rémillard (2001), A nonparametric test
  of serial independence for times series and residuals. \emph{Journal
    of Multivariate Analysis}, 79:191--218.
  
  I. Kojadinovic and J. Yan (2011), Tests of multivariate serial
  independence based on a Möbius decomposition of the independence
  empirical copula process, \emph{Annals of the Institute of Statistical
  Mathematics}, 63:347--373. 
}

\seealso{
  \code{\link{serialIndepTest}},
  \code{\link{indepTest}},
  \code{\link{multIndepTest}},
  \code{\link{dependogram}}
}
\examples{
## A multivariate time series
d <- 2
n <- 100
param <- 0.25
ar <- matrix(0,2*n,d)
ar[1,] <- rnorm(d)
for (i in 2:(2*n))
  ar[i,] <- matrix(param,d,d) \%*\% ar[i-1,] + rnorm(d)
x <- ar[(n+1):(2*n),]

## Run the test
test <- multSerialIndepTest(x,3)
test

## Display the dependogram
dependogram(test,print=TRUE)
}

\encoding{UTF-8}
\keyword{htest}
