\name{gofCopula}
\alias{gofCopula}

\title{Goodness-of-fit tests for copulas}
\description{
  Goodness-of-fit tests for copulas based on the empirical process
  comparing the so-called \emph{empirical copula} with a
  parametric estimate of the copula derived under the null hypothesis.
  Approximate p-values for the test statistic can be obtained either
  using the so-called \emph{parametric bootstrap} (see the two first
  references) or by means of a fast \emph{multiplier} approach (see the two
  last references).  
}
\usage{
gofCopula(copula, x, N = 10000, estimation = "likelihood",
          method = "parametric.bootstrap", M = 2500)
}

\arguments{
  \item{copula}{ object of class \code{"copula"} representing the
    hypothesized copula family.}
  \item{x}{ a data matrix. }
  \item{N}{ number of multiplier or bootstrap iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{estimation}{ estimation method to be used to estimate the
    dependence parameter(s); can be either \code{"likelihood"},
    \code{"kendall"} or \code{"spearman"}. }
  \item{method}{ method for generating realizations of the test
    statistic under the null hypothesis; can be either
    \code{"multiplier"} or \code{"parametric.bootstrap"}. }
  \item{M}{ size of the Monte Carlo integration sample; see the two last
    references. }
}
\details{
  Thus far, the multiplier approach is implemented for six copula
  families: the Clayton, Gumbel, Frank, Plackett, normal and t. For the
  normal and t copulas, several dependence structures can be
  hypothesized: \code{"ex"} for exchangeable, \code{"ar1"} for AR(1),
  \code{"toep"} for  Toeplitz, and \code{"un"} for unstructured (see
  \code{\link{ellipCopula}}). For the t copula, \code{df.fixed} has to
  be set to \code{TRUE}.
  
  In the bivariate one-parameter case, the dependence parameter of
  the hypothesized copula family can be estimated either by maximizing
  the pseudo-likelihood or by inverting Kendall's tau or Spearman's rho. In
  the multivariate (multiparameter) case, only maximum pseudo-likelihood
  estimation can be used.
  
  If \code{"estimation"} is set to \code{"likelihood"}, the parameter values
  of \code{"copula"} are used as starting values for the pseudo-likelihood
  maximization.}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistics. }
  \item{pvalue}{ corresponding approximate p-value. }
}
\references{

  C. Genest and B. Remillard (2008). Validity of the parametric
  bootstrap for goodness-of-fit testing in semiparametric models.
  \emph{Annales de l'Institut Henri Poincare: Probabilites et Statistiques},
  44, in press.

  C. Genest, B. Remillard and D. Beaudoin (2008). Goodness-of-fit tests
  for copulas: A review and a power study. \emph{Insurance: Mathematics and
    Economics}, 42, in press.
  
  I. Kojadinovic and J. Yan (2008). Fast goodness-of-fit tests for
  copulas. Submitted.

  I. Kojadinovic and J. Yan (2008). A goodness-of-fit test for
  multivariate multiparameter copulas based on multiplier central limit
  theorems. Submitted.
}
\seealso{ \code{\link{fitCopula}}, \code{\link{ellipCopula}}. }
\examples{
## A two-dimensional data example
x <- rcopula(claytonCopula(3), 200)

\dontrun{
## Taking too long with method = "parametric.bootstrap"
## Does the Gumbel family seem to be a good choice? 
gofCopula(gumbelCopula(1), x, estimation="kendall")
## What about the Clayton family? 
gofCopula(claytonCopula(1), x, estimation="kendall")
## Pseudo-likelihood estimation where starting value are obtained
## by inverting Kendall's tau
clayton.cop <- claytonCopula(1)
clayton.cop@parameters[1] <- calibKendallsTau(clayton.cop, cor(x, method="kendall")[1,2])
gofCopula(clayton.cop, x)
}
## A three-dimensional example
x <- rcopula(normalCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"),200)

\dontrun{
## Taking too long with method = "parametric.bootstrap"
## Does the Clayton family seem to be a good choice? 
gofCopula(gumbelCopula(1, dim = 3), x)
## What about the normal copula?
normal.copula <- normalCopula(rep(0, 3), dim = 3, dispstr = "un")
gofCopula(normal.copula, x)
}
}

\keyword{htest}
\keyword{models}
\keyword{multivariate}
