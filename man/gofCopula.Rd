\name{gofCopula}
\alias{gofCopula}
\title{Goodness-of-fit Tests for Copulas}
\description{
  Goodness-of-fit tests for copulas based on the empirical process
  comparing the empirical copula with a
  parametric estimate of the copula derived under the null hypothesis.
  The test statistic is the Cramer-von Mises functional
  \eqn{\mathrm{S_n}}{S[n]} defined
  in Equation (2) of Genest, Remillard and Beaudoin (2009).
  Approximate p-values for the test statistic can be obtained either
  using the \emph{parametric bootstrap} (see the two first
  references) or by means of a fast \emph{multiplier} approach (see the
  last two references).
}
\usage{
gofCopula(copula, x, N = 1000,
          method = c("Sn", "AnChisq", "AnGamma", "SnB", "SnC"),
          estim.method = c("mpl", "ml", "irho", "itau"),
          simulation = c("pb", "mult"), verbose = TRUE, print.every = NULL,
          optim.method = "BFGS", optim.control = list(maxit=20))
}
\arguments{
  \item{copula}{ object of class \code{"\linkS4class{copula}"} representing the
    hypothesized copula family.}
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{ number of bootstrap or multiplier iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{method}{a \code{\link{character}} string specifying the
    goodness-of-fit test statistic to be used, which has to be one (or a unique
    abbreviation) of \code{"AnChisq"}, \code{"AnGamma"}, \code{"SnB"},
    \code{"SnC"}, or \code{"Sn"}.

    \code{"Sn"} is \eqn{\mathrm{S_n}}{S[n]} from Equation (2) in Genest,
    R\enc{é}{e}millard, Beaudoin (2009). For the others see
    \code{method} in \code{\link{gofTstat}}.
  }
  \item{estim.method}{a character string specifying the estimation method to
    be used to estimate the dependence parameter(s); can be either \code{"mpl"}
    (maximum pseudo-likelihood), \code{"itau"} (inversion of
    Kendall's tau) or \code{"irho"} (inversion of Spearman's rho).}
  \item{simulation}{a string specifying the simulation method for
    generating realizations of the test statistic under the null
    hypothesis; can be either \code{"pb"} (parametric bootstrap) or
    \code{"mult"} (multiplier).}
  \item{print.every}{ is deprecated in favor of \code{verbose}.}
  \item{verbose}{a logical specifying if progress of the bootstrap
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{optim.method, optim.control}{the \code{method} and
    \code{control} arguments for \code{\link{optim}()}, see there.}
}
\details{
  If the parametric bootstrap is used, the dependence parameters of
  the hypothesized copula family can be estimated either by maximizing
  the pseudo-likelihood, by inverting Kendall's tau, or by inverting
  Spearman's rho.  If the multiplier is used, any estimation method
  can be used in the bivariate case, but only maximum pseudo-likelihood
  estimation can be used in the multivariate (multiparameter) case.

  For the normal and t copulas, several dependence structures can be
  hypothesized: \code{"ex"} for exchangeable, \code{"ar1"} for AR(1),
  \code{"toep"} for Toeplitz, and \code{"un"} for unstructured (see
  \code{\link{ellipCopula}}). For the t copula, \code{"df.fixed"} has to
  be set to \code{TRUE}, which implies that the degrees of freedom are
  not considered as a parameter to be estimated.

  Thus far, the multiplier approach is implemented for six copula
  families: the Clayton, Gumbel, Frank, Plackett, normal and t.

  Although the processes involved in the multiplier and the parametric
  bootstrap-based test are asymptotically equivalent under the null,
  note that the finite-sample behavior of the two tests might differ
  significantly.

  Also note that in case of the parametric bootstrap, a finite-sample
  correction of type \eqn{(sum(T0 >= T)+0.5)/(N+1)} is made, where
  \eqn{T} and \eqn{T0} are the test statistic and bootstrapped test
  statistc, respectively.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
  \item{parameter}{ estimates of the parameters for the hypothesized
    copula family. }
}
\references{
  Genest, C. and R\enc{é}{e}millard, B. (2008). Validity of the parametric
  bootstrap for goodness-of-fit testing in semiparametric models.
  \emph{Annales de l'Institut Henri Poincare: Probabilites et Statistiques}
  \bold{44}, 1096--1127.

  Genest, C., R\enc{é}{e}millard, B., and Beaudoin,
  D. (2009). Goodness-of-fit tests for copulas: A review and a power study.
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--214.

  Kojadinovic, I., Yan, J., and Holmes M. (2011).
  Fast large-sample goodness-of-fit tests for copulas.
  \emph{Statistica Sinica} \bold{21}, 841--871.

  Kojadinovic, I. and Yan, J. (2011). A goodness-of-fit test for
  multivariate multiparameter copulas based on multiplier central limit
  theorems. \emph{Statistics and Computing} \bold{21}, 17--30.

  Kojadinovic, I. and Yan, J. (2010).
  Modeling Multivariate Distributions with Continuous Margins Using the
  copula R Package.
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
  \code{\link{fitCopula}()}, \code{\link{ellipCopula}}.

  \code{\link{gnacopula}} for other goodness-of-fit tests for (nested)
  Archimedean copulas.
}
\examples{
## the following example is available in batch through
## demo(gofCopula)
%% Note that  ../tests/gof-ex.R (at the end) has "proof-of-concept" tests
\dontrun{
## A two-dimensional data example ----------------------------------
x <- rCopula(200, claytonCopula(3))

tau. <- cor(x, method="kendall")
## Does the Gumbel family seem to be a good choice?
thG <- iTau(gumbelCopula(), tau.)
gofCopula(gumbelCopula(thG), x)
# SnC: really s..l..o..w.. --- SnB is *EVEN* slower
gofCopula(gumbelCopula(thG), x, method = "SnC")
## What about the Clayton family?
thC <- iTau(claytonCopula(), tau.)
gofCopula(claytonCopula(thC), x)
gofCopula(claytonCopula(thC), x, method = "AnChisq")

## The same with a different estimation method
gofCopula(gumbelCopula (thG), x, estim.method="itau")
gofCopula(claytonCopula(thC), x, estim.method="itau")



## A three-dimensional example  ------------------------------------
x <- rCopula(200, tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"))

## Does the Clayton family seem to be a good choice?
gofCopula(gumbelCopula(1, dim = 3), x)
## What about the t copula?
t.copula <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)
## this is *VERY* slow currently %% FIXME ??
gofCopula(t.copula, x)

## The same with a different estimation method
gofCopula(gumbelCopula(1, dim = 3), x, estim.method="itau")
gofCopula(t.copula,                 x, estim.method="itau")

## The same using the multiplier approach
gofCopula(gumbelCopula(1, dim = 3), x, simulation="mult")
gofCopula(t.copula,                 x, simulation="mult")
}% dontrun
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
