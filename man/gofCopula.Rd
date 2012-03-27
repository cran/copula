\name{gofCopula}
\alias{gofCopula}

\title{Goodness-of-fit tests for copulas}
\description{
  Goodness-of-fit tests for copulas based on the empirical process
  comparing the empirical copula with a
  parametric estimate of the copula derived under the null hypothesis.
  The test statistic is the Cramer-von Mises functional Sn defined
  in equation (2) of Genest, Remillard and Beaudoin (2009).
  Approximate p-values for the test statistic can be obtained either
  using the \emph{parametric bootstrap} (see the two first
  references) or by means of a fast \emph{multiplier} approach (see the two
  last references).
}
\usage{
gofCopula(copula, x, N = 1000, method = "mpl",
          simulation = c("pb", "mult"),
          print.every = 100, optim.method = "BFGS")
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
  \item{print.every}{ progress is printed every \code{"print.every"}
    iterations. No progress is printed if it is nonpositive.}
  \item{optim.method}{ the method for \code{\link{optim}()}.}
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

  %Although the processes involved in the multiplier and the parametric
  %bootstrap-based test are asymptotically equivalent, the finite-sample
  %behavior of the two tests might differ significantly.
}
\value{
  An object of \code{\link{class}} \code{"gofCopula"} which is a list
  with components
  \item{statistic}{ value of the test statistic. }
  \item{pvalue}{ corresponding approximate p-value. }
  \item{parameters}{ estimates of the parameters for the hypothesized
    copula family. }
}
\references{
  C. Genest and B. Remillard (2008). Validity of the parametric
  bootstrap for goodness-of-fit testing in semiparametric models.
  \emph{Annales de l'Institut Henri Poincare: Probabilites et Statistiques}
  \bold{44}, 1096--1127.

  C. Genest, B. Remillard and D. Beaudoin (2009). Goodness-of-fit tests
  for copulas: A review and a power study.
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--214.

  I. Kojadinovic, J. Yan and M. Holmes (2011).
  Fast large-sample goodness-of-fit tests for copulas.
  \emph{Statistica Sinica} \bold{21}, 841--871.

  I. Kojadinovic and J. Yan (2011). A goodness-of-fit test for
  multivariate multiparameter copulas based on multiplier central limit
  theorems. \emph{Statistics and Computing} \bold{21}, 17--30.

  I. Kojadinovic and J. Yan (2010).
  Modeling Multivariate Distributions with Continuous Margins Using the copula R Package.
  \emph{Journal of Statistical Software} \bold{34}(9), 1--20.
  \url{http://www.jstatsoft.org/v34/i09/}.
}
\note{
  These tests were derived under the assumption of continuous margins,
  which implies that ties occur with probability zero. The
  presence of ties in the data might substantially affect the
  approximate p-values. One way of dealing with ties was suggested in the
  last reference.
}
\seealso{
  \code{\link{fitCopula}}, \code{\link{ellipCopula}}.

  \code{\link{gnacopula}} for other goodness-of-fit tests for (nested)
  Archimedean copulas.
}
\examples{
## the following example is available in batch through
## demo(gofCopula)
%% FIXME: We need some *run* examples, either here or in ../tests/ !!
\dontrun{
## A two-dimensional data example ----------------------------------
x <- rcopula(claytonCopula(3), 200)

## Does the Gumbel family seem to be a good choice?
gofCopula(gumbelCopula(1), x)
## What about the Clayton family?
gofCopula(claytonCopula(1), x)

## The same with a different estimation method
gofCopula(gumbelCopula (1), x, method="itau")
gofCopula(claytonCopula(1), x, method="itau")


## A three-dimensional example  ------------------------------------
x <- rcopula(tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"),
             200)

## Does the Clayton family seem to be a good choice?
gofCopula(gumbelCopula(1, dim = 3), x)
## What about the t copula?
t.copula <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)
## this is *VERY* slow currently %% FIXME ??
gofCopula(t.copula, x)

## The same with a different estimation method
gofCopula(gumbelCopula(1, dim = 3), x, method="itau")
gofCopula(t.copula,                 x, method="itau")

## The same using the multiplier approach
gofCopula(gumbelCopula(1, dim = 3), x, simulation="mult")
gofCopula(t.copula,                 x, simulation="mult")
}% dontrun

}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
