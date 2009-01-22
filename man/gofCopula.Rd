\name{gofCopula}
\alias{gofCopula}

\title{Goodness-of-fit tests for copulas}
\description{
  Goodness-of-fit tests for copulas based on the empirical process
  comparing the \emph{empirical copula} with a
  parametric estimate of the copula derived under the null hypothesis.
  Approximate p-values for the test statistic can be obtained either
  using the \emph{parametric bootstrap} (see the two first
  references) or by means of a fast \emph{multiplier} approach that can
  be used when the parametric bootstrap is too slow (see the two
  last references).  
}
\usage{
gofCopula(copula, x, N = 1000, method = "mpl",
          simulation = "pb", grid = "h0", R = 10,
          m = nrow(x), G = nrow(x), M = 2500,
          print.every = 100, optim.method = "Nelder-Mead")
}

\arguments{
  \item{copula}{ object of class \code{"copula"} representing the
    hypothesized copula family.}
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{ number of bootstrap or multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{method}{ estimation method to be used to estimate the
    dependence parameter(s); can be either \code{"mpl"}
    (maximum pseudo-likelihood), \code{"itau"} (inversion of
    Kendall's tau) or \code{"irho"} (inversion of Spearman's rho). }
  \item{simulation}{ simulation method for generating realizations
    of the test statistic under the null hypothesis; can be either
    \code{"pb"} (parametric bootstrap) or  \code{"mult"} (multiplier).}
  \item{grid}{ for simulation method \code{"mult"}, grid points at which the
    goodness-of-fit process is evaluated; can be either \code{"h0"} (the
    goodness-of-fit process is evaluated at points randomly generated
    from the hypothesized copula) or \code{"po"} (the
    goodness-of-fit process is evaluated at the available
    pseudo-observations); see the two last references. }
  \item{R}{ for simulation method \code{"mult"}, number of replications of the
    basic test; should be set to 1 for very large samples. }
  \item{m}{ for simulation method \code{"mult"}, size of the sample
    used to compute the influence functions; see the two last references. }
  \item{G}{ for simulation method \code{"mult"}, size of the grid
    if \code{"grid"} is set to \code{"h0"}; see the two last references. }
  \item{M}{ for simulation method \code{"mult"}, size of the Monte Carlo
    integration sample; see the two last references. }
  \item{print.every}{ progress is printed every \code{"print.every"}
    iterations. }
  \item{optim.method}{ the method for \code{"\link{optim}"}. }
}
\details{
  
  If the parametric bootstrap is used, the dependence parameters of
  the hypothesized copula family can be estimated either by maximizing
  the pseudo-likelihood or by inverting Kendall's tau or Spearman's rho.
  If the multiplier is used, any estimation method can be used in the
  bivariate case, but only maximum pseudo-likelihood
  estimation can be used in the multivariate (multiparameter) case.
  
  For the normal and t copulas, several dependence structures can be
  hypothesized: \code{"ex"} for exchangeable, \code{"ar1"} for AR(1),
  \code{"toep"} for  Toeplitz, and \code{"un"} for unstructured (see
  \code{\link{ellipCopula}}). For the t copula, \code{"df.fixed"} has to
  be set to \code{TRUE}, which implies that the degrees of freedom are
  not considered as a parameter to be estimated.

  Thus far, the multiplier approach is implemented for six copula
  families: the Clayton, Gumbel, Frank, Plackett, normal and t. 
  
  The parameter \code{"R"}, used when the simulation method is \code{"mult"},
  specifies the number of times that the multiplier test will be
  repeated. Indeed, when the available sample is smalle, the multiplier
  approach can show poor repeatability with respect to the returned p-value (see the last two
  references). As the sample size increases, the value of
  \code{"R"} can be decreased.

  The parameter \code{"grid"}, used when the simulation method is
  \code{"mult"}, should be set as follows: if the sample size is smaller
  than 500, \code{"grid"} should be set to \code{"h0"} (otherwise the
  test will be too liberal); when the sample size is greater than 500,
  \code{"grid"} should be set to \code{"po"}, as this decreases
  repeatability problems. 
  
  Although the processes involved in the multiplier and the parametric
  bootstrap-based test are asymptotically equivalent, the finite-sample
  behavior of the two tests might differ significantly.
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic, or median of the test
    statistics if \code{"R"} multiplier replications are performed. }
  \item{pvalue}{ corresponding approximate p-value, or median of the
    p-values if \code{"R"} multiplier replications are performed. }
  \item{sd.pvalue}{ standard deviation of the pvalues when
   the multiplier test is replicated \code{"R"} times. }
  \item{parameters}{ estimates of the parameters for the hypothesized
    copula family. }
}
\references{

  C. Genest and B. Remillard (2008). Validity of the parametric
  bootstrap for goodness-of-fit testing in semiparametric models.
  \emph{Annales de l'Institut Henri Poincare: Probabilites et Statistiques},
  44, 1096-1127.

  C. Genest, B. Remillard and D. Beaudoin (2008). Goodness-of-fit tests
  for copulas: A review and a power study. \emph{Insurance: Mathematics and
    Economics}, 44, in press.
  
  I. Kojadinovic and J. Yan (2008). Fast large-sample goodness-of-fit
  tests for copulas. Submitted.

  I. Kojadinovic and J. Yan (2008). A goodness-of-fit test for
  multivariate multiparameter copulas based on multiplier central limit
  theorems. Submitted.
}
\seealso{ \code{\link{fitCopula}}, \code{\link{ellipCopula}}. }
\examples{
## the following example is available in batch through
## demo(gofCopula)
\dontrun{
## A two-dimensional data example
x <- rcopula(claytonCopula(3), 200)

## Does the Gumbel family seem to be a good choice? 
gofCopula(gumbelCopula(1), x) 
## What about the Clayton family? 
gofCopula(claytonCopula(1), x)

## The same with a different estimation method
gofCopula(gumbelCopula(1), x, method="itau")
gofCopula(claytonCopula(1), x, method="itau")


## A three-dimensional example
x <- rcopula(tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"),200)

## Does the Clayton family seem to be a good choice? 
gofCopula(gumbelCopula(1, dim = 3), x)
## What about the t copula?
t.copula <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)
gofCopula(t.copula, x) 

## The same with a different estimation method
gofCopula(gumbelCopula(1, dim = 3), x, method="itau")
gofCopula(t.copula, x, method="itau")

## The same using the multiplier approach
gofCopula(gumbelCopula(1, dim = 3), x, simulation="mult")
gofCopula(t.copula, x, simulation="mult")
}
}

\keyword{htest}
\keyword{models}
\keyword{multivariate}
