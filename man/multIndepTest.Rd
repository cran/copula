\name{multIndepTest}
\alias{multIndepTest}

\title{Independence test among continuous random vectors
  based on the empirical copula process}

\description{
  Analog of the independence test based on the empirical copula process
  proposed by Christian Genest and Bruno Rémillard (see
  \code{\link{indepTest}}) for \emph{random vectors}. The main
  difference comes from the fact that critical values and p-values are
  obtainted through the bootstrap/permutation methodology, since, here,
  test statistics are not distribution-free.
}

\usage{
multIndepTest(x, d, m=length(d), N=1000, alpha=0.05, print.every = 100)
}

\arguments{
  \item{x}{Data frame or data matrix containing realizations (one per
    line) of the random vectors whose independence is to be tested.}
  \item{d}{Dimensions of the random vectors whose realizations are given
    in \code{x}. It is required that \code{sum(d)=ncol(x)}.}
  \item{m}{Maximum cardinality of the subsets of random vectors for
    which a test statistic is to be computed. It makes sense to consider
    \code{m << p} especially when \code{p} is large.}
  \item{N}{Number of bootstrap/permutation samples.}
  \item{alpha}{Significance level used in the computation of the
    critical values for the test statistics.}
  \item{print.every}{Progress is printed every \code{"print.every"}
    iterations. No progress is printed if it is nonpositive.}
}
\details{
  See the references below for more details, especially the last one.
}

\value{
  The function \code{"multIndepTest"} returns an object of class
  \code{"indepTest"} whose attributes are: \code{subsets},
  \code{statistics}, \code{critical.values}, \code{pvalues},
  \code{fisher.pvalue} (a p-value resulting from a combination \emph{à la}
  Fisher of the subset statistic p-values), \code{tippett.pvalue} (a p-value
  resulting from a combination \emph{à la} Tippett of the subset
  statistic p-values), \code{alpha} (global significance level of the
  test), \code{beta} (1 - \code{beta} is the significance level per
  statistic), \code{global.statistic} (value of the global Cramér-von
  Mises statistic derived directly from the independence empirical
  copula process - see \code{In} in the last reference) and
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
  
  I. Kojadinovic and M. Holmes (2009), Tests of independence among
  continuous random vectors based on Cramér-von Mises functionals
  of the empirical copula process. \emph{Journal of Multivariate
    Analysis}, 100:1137--1154.
}

\seealso{
  \code{\link{indepTest}},
  \code{\link{serialIndepTest}},
  \code{\link{multSerialIndepTest}},
  \code{\link{dependogram}}.}

\examples{
## Consider the following example taken from
## Kojadinovic and Holmes (2008):

n <- 100

## Generate data
y <- matrix(rnorm(6*n),n,6)
y[,1] <- y[,2]/2 + sqrt(3)/2*y[,1]
y[,3] <- y[,4]/2 + sqrt(3)/2*y[,3]
y[,5] <- y[,6]/2 + sqrt(3)/2*y[,5]
    
nc <- normalCopula(0.3,dim=3)
x <- cbind(y,rcopula(nc,n),rcopula(nc,n))
       
x[,1] <- abs(x[,1]) * sign(x[,3] * x[,5])
x[,2] <- abs(x[,2]) * sign(x[,3] * x[,5])
x[,7] <- x[,7] + x[,10]
x[,8] <- x[,8] + x[,11]
x[,9] <- x[,9] + x[,12]

## Dimensions of the random vectors
d <- c(2,2,2,3,3)

## Run the test
test <- multIndepTest(x,d)
test

## Display the dependogram
dependogram(test,print=TRUE)
}

\encoding{UTF-8}
\keyword{htest}
