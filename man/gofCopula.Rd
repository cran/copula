\name{gofCopula}
\title{Goodness-of-fit Tests for Copulas}
\alias{gofCopula}
\alias{gofCopula-methods}
\alias{gofCopula,copula-method}
\alias{gofCopula,parCopula-method}
\alias{gofCopula,rotCopula-method}
\alias{gofMB}
\alias{gofPB}
\description{
  The goodness-of-fit tests are based, by default, on the empirical
  process comparing the empirical copula with a parametric estimate of
  the copula derived under the null hypothesis, the default test
  statistic, "Sn", being the Cramer-von Mises functional \eqn{S_n}{S[n]}
  defined in Equation (2) of Genest, Remillard and Beaudoin (2009). In
  that case, approximate p-values for the test statistic can be obtained
  either using a \emph{parametric bootstrap} (see references two and
  three) or by means of a faster \emph{multiplier} approach (see
  references four and five).

  Alternative test statistics can be used, in particular if a
  \emph{parametric bootstrap} is employed.

  The prinicipal function is \code{gofCopula()} which, depending on
  \code{simulation} either calls \code{gofPB()} or \code{gofMB()}.
}
\usage{
## Generic [and "rotCopula" method] ------ Main function ------
gofCopula(copula, x, \dots)
\S4method{gofCopula}{copula}(copula, x, N = 1000,
          method = c("Sn", "SnB", "SnC", "Rn"),
          estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
          simulation = c("pb", "mult"), test.method = c("family", "single"),
          verbose = interactive(), ties = NA,
          ties.method = c("max", "average", "first", "last", "random", "min"),
          fit.ties.meth = eval(formals(rank)$ties.method), \dots)

## (Deprecated) internal 'helper' functions : ---
gofPB(copula, x, N, method = c("Sn", "SnB", "SnC"),
      estim.method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
      trafo.method = if(method == "Sn") "none" else c("cCopula", "htrafo"),
      trafoArgs = list(), test.method = c("family", "single"),
      verbose = interactive(), useR = FALSE, ties = NA,
      ties.method = c("max", "average", "first", "last", "random", "min"),
      fit.ties.meth = eval(formals(rank)$ties.method), \dots)

gofMB(copula, x, N, method = c("Sn", "Rn"),
      estim.method = c("mpl", "ml", "itau", "irho"),
      test.method = c("family", "single"), verbose = interactive(),
      useR = FALSE, m = 1/2, zeta.m = 0, b = 1/sqrt(nrow(x)),
      ties.method = c("max", "average", "first", "last", "random", "min"),
      fit.ties.meth = eval(formals(rank)$ties.method), \dots)
} % NB: keep internal helper functions on this help page, even if hidden !
\arguments{
  \item{copula}{object of class \code{"\linkS4class{copula}"} representing the
    hypothesized copula family.}
  \item{x}{a data matrix that will be transformed to pseudo-observations
    using \code{\link{pobs}()}.}
  \item{N}{number of bootstrap or multiplier replications to be used to
    obtain approximate realizations of the test statistic under the null
    hypothesis.}
  \item{method}{a \code{\link{character}} string specifying the
    goodness-of-fit test statistic to be used.  For \code{simulation = "pb"},
    one of "Sn", "SnB" or "SnC" with \code{trafo.method != "none"} if
      \code{method != "Sn"}.
    For \code{simulation = "mult"}, one of
    \code{"Sn"} or \code{"Rn"}, where the latter is \eqn{R_n} from
    Genest et al. (2013).}
  \item{estim.method}{a \code{\link{character}} string specifying the
    estimation method to be used to estimate the dependence
    parameter(s); see \code{\link{fitCopula}()}.}
  \item{simulation}{a string specifying the resampling method for
    generating approximate realizations of the test statistic under the null
    hypothesis; can be either \code{"pb"} (parametric bootstrap) or
    \code{"mult"} (multiplier).}
  \item{test.method}{a \code{\link{character}} string specifying the
    test method to be used. Only in exceptional cases should this be
    different from the default \code{test.method = "family"}.
    If \code{test.method = "single"}, a test precisely for the provided
    copula (not its parametric family) is conducted. This makes sense only for
    specific applications such as testing random number generators.}
  \item{verbose}{a logical specifying if progress of the parametric bootstrap
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{\dots}{for \code{gofCopula}, additional arguments passed to
    \code{gofPB()} or \code{gofMB()};
    for \code{gofPB()} and \code{gofMB()}: additional arguments passed
    to \code{\link{fitCopula}()}.  These may notably contain
    \code{hideWarnings}, and
    \code{optim.method}, \code{optim.control}, \code{lower},
    or \code{upper} depending on the \code{optim.method}.}
  \item{trafo.method}{only for the parametric bootstrap (\code{"pb"}):
    String specifying the transformation to
    \eqn{U[0,1]^d}; either \code{"none"} or one of \code{"cCopula"},
    see \code{\link{cCopula}()}, or \code{"htrafo"}, see
    \code{\link{htrafo}()}.  If \code{method != "Sn"}, one needs to set
    \code{trafo.method != "none"}.}
  \item{trafoArgs}{only for the parametric bootstrap.  A
    \code{\link{list}} of optional arguments passed to the
    transformation method (see \code{trafo.method} above).}
  \item{useR}{logical indicating whether an \R or C implementation is
    used.}
  \item{ties.method}{string specifying how ranks should be computed,
      except for fitting, if there are ties in any of the coordinate
      samples of \code{x}; passed to \code{\link{pobs}}.}
  \item{fit.ties.meth}{string specifying how ranks should be computed
    when fitting by maximum pseudo-likelihood if there are ties in any
    of the coordinate samples of \code{x}; passed to \code{\link{pobs}}.}
  \item{ties}{only for the parametric bootstrap.  Logical indicating
    whether a version of the parametric boostrap adapted to the
    presence of ties in any of the coordinate samples of \code{x}
    should be used; the default value of \code{NA} indicates that the
    presence/absence of ties will be checked for automatically.}
  \item{m, zeta.m}{only for the multiplier with \code{method = "Rn"}.
    \code{m} is the power and \code{zeta.m} is the adjustment
    parameter \eqn{\zeta_m}{zeta.m} for the denominator of the test
    statistic.}
  \item{b}{only for the multiplier.  \code{b} is the bandwidth required
    for the estimation of the first-order partial derivatives based on
    the empirical copula.}
}
\details{
  If the parametric bootstrap is used, the dependence parameters of
  the hypothesized copula family can be estimated by any estimation
  method available for the family, up to a few exceptions.
  If the multiplier is used, any of the rank-based methods
  can be used in the bivariate case, but only maximum pseudo-likelihood
  estimation can be used in the multivariate (multiparameter) case.

  The price to pay for the higher computational efficiency of the
  multiplier is more programming work as certain
  partial derivatives need to be computed for each hypothesized
  parametric copula family. When estimation is based on maximization of
  the pseudo-likelihood, these have been implemented for six copula
  families thus far: the Clayton, Gumbel-Hougaard, Frank, Plackett,
  normal and \eqn{t} copula families.  For other families, numerical
  differentiation based on \code{\link[numDeriv]{grad}()} from package
  \CRANpkg{numDeriv} is used (and a warning message is displayed).

  Although the empirical processes involved in the multiplier and the
  parametric bootstrap-based test are asymptotically equivalent under
  the null, the finite-sample behavior of the two tests might differ
  significantly.

  Both for the parametric bootstrap and the multiplier,
  the approximate p-value is computed as
  \deqn{(0.5 +\sum_{b=1}^N\mathbf{1}_{\{T_b\ge T\}})/(N+1),}{%
        (0.5 + sum(T[b] >= T, b=1, .., N)) / (N+1),}
  where \eqn{T} and \eqn{T_b}{T[b]} denote the test statistic and
  the bootstrapped test statistc, respectively. This ensures that the
  approximate p-value is a number strictly between 0 and 1, which is
  sometimes necessary for further treatments. See Pesarin (2001) for
  more details.

  For the normal and \eqn{t} copulas, several dependence structures can
  be hypothesized: \code{"ex"} for exchangeable, \code{"ar1"} for AR(1),
  \code{"toep"} for Toeplitz, and \code{"un"} for unstructured (see
  \code{\link{ellipCopula}()}). For the \eqn{t} copula,
  \code{"df.fixed"} has to be set to \code{TRUE}, which implies that the
  degrees of freedom are not considered as a parameter to be estimated.

  The former argument \code{print.every} is deprecated and not
  supported anymore; use \code{verbose} instead.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{value of the test statistic.}
  \item{p.value}{corresponding approximate p-value.}
  \item{parameter}{estimates of the parameters for the hypothesized
    copula family.}
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Genest, C., Huang, W., and Dufour, J.-M. (2013).  A regularized
  goodness-of-fit test for copulas.  \emph{Journal de la
  Soci\enc{é}{e}t\enc{é}{e} fran\enc{ç}{c}aise de statistique}
  \bold{154}, 64--77.

  Genest, C. and R\enc{é}{e}millard, B. (2008). Validity of the parametric
  bootstrap for goodness-of-fit testing in semiparametric models.
  \emph{Annales de l'Institut Henri Poincare: Probabilites et Statistiques}
  \bold{44}, 1096--1127.

  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009).
  Goodness-of-fit tests for copulas: A review and a power study.
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
  \emph{Journal of Statistical Software} \bold{34}(9), 1--20,
  \url{http://www.jstatsoft.org/v34/i09/}.

  Kojadinovic, I. (2017). Some copula inference procedures adapted to
  the presence of ties. \emph{Computational Statistics and Data
  Analysis} \bold{112}, 24--41, \url{http://arxiv.org/abs/1609.05519}.

  Pesarin, F. (2001).
  \emph{Multivariate Permutation Tests: With Applications in Biostatistics}.
  Wiley.
}
\note{
  These tests were theoretically studied and implemented under the
  assumption of continuous margins, which implies that ties in the
  component samples occur with probability zero. The presence of ties in
  the data might substantially affect the approximate p-values. Through
  argument \code{ties}, the user can however select a version of the
  parametric bootstrap adapted to the presence of ties. No such adaption
  exists for the multiplier for the moment.

  %% Since \R is widely used by practitioners, a word of warning concerning
  %% goodness-of-fit tests \emph{in general} is also advisable.
  %% Goodness-of-fit tests are often (ab)used in practice to
  %% \dQuote{justify} an assumption under which one then continues to work
  %% (carelessly).  From a mathematical point of view, this is not correct.
}
\seealso{
  \code{\link{fitCopula}()} for the underlying estimation procedure and
  \code{\link{gofTstat}()} for details on *some* of the available test
  statistics.
}

%% Examples == ../demo/gofCopula.R __keep >> EXACTLY << in sync!__
%% Note that  ../tests/gof-ex.R (at the end) has "proof-of-concept"
%% tests

\examples{
## The following example is available in batch through
## demo(gofCopula)
\donttest{
n <- 200; N <- 1000 # realistic (but too large for interactive use)
n <-  60; N <-  200 # (time (and tree !) saving ...)

## A two-dimensional data example ----------------------------------
set.seed(271)
x <- rCopula(n, claytonCopula(3))

## Does the Gumbel family seem to be a good choice (statistic "Sn")?
gofCopula(gumbelCopula(), x, N=N)
## With "SnC", really s..l..o..w.. --- with "SnB", *EVEN* slower
gofCopula(gumbelCopula(), x, N=N, method = "SnC", trafo.method = "cCopula")
## What about the Clayton family?
gofCopula(claytonCopula(), x, N=N)

## Similar with a different estimation method
gofCopula(gumbelCopula (), x, N=N, estim.method="itau")
gofCopula(claytonCopula(), x, N=N, estim.method="itau")


## A three-dimensional example  ------------------------------------
x <- rCopula(n, tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"))

## Does the Gumbel family seem to be a good choice?
g.copula <- gumbelCopula(dim = 3)
gofCopula(g.copula, x, N=N)
## What about the t copula?
t.copula <- tCopula(dim = 3, dispstr = "un", df.fixed = TRUE)
if(FALSE) ## this is *VERY* slow currently
  gofCopula(t.copula, x, N=N)

## The same with a different estimation method
gofCopula(g.copula, x, N=N, estim.method="itau")
if(FALSE) # still really slow
  gofCopula(t.copula, x, N=N, estim.method="itau")

## The same using the multiplier approach
gofCopula(g.copula, x, N=N, simulation="mult")
gofCopula(t.copula, x, N=N, simulation="mult")
if(FALSE) # no yet possible
    gofCopula(t.copula, x, N=N, simulation="mult", estim.method="itau")
}% dont..
}
\keyword{htest}
\keyword{models}
\keyword{goodness-of-fit}
\keyword{multivariate}
