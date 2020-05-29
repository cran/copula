\name{fitCopula}
\title{Fitting Copulas to Data -- Copula Parameter Estimation}
\alias{fitCopula}
\alias{fitCopula-methods}
\alias{fitCopula,copula-method}
\alias{fitCopula,parCopula-method}
\alias{fitCopula,rotCopula-method}
\alias{loglikCopula}
\alias{optimMeth}
\description{
  Parameter estimation of copulas, i.e., fitting of a copula model to
  multivariate (possibly \dQuote{pseudo}) observations.
}
\usage{
loglikCopula(param = getTheta(copula), u, copula,
             error = c("-Inf", "warn-Inf", "let-it-be"))

## Generic [and "rotCopula" method] :
fitCopula(copula, data, \dots)
\S4method{fitCopula}{copula}(copula, data,
          method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
          posDef = is(copula, "ellipCopula"),
          start = NULL, lower = NULL, upper = NULL,
          optim.method = optimMeth(copula, method, dim = d),
          optim.control = list(maxit=1000),
          estimate.variance = NA, hideWarnings = FALSE, \dots)

optimMeth(copula, method, dim)
}
\arguments{
  \item{param}{vector of \emph{free} (see \code{\link{isFree}()} and
    \code{\link{getTheta}()}) parameter values.}
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations in
    \eqn{[0,1]^d} for computing the copula log-likelihood, where \eqn{n} denotes
    the sample size and \eqn{d} the dimension.  Consider applying the function
    \code{\link{pobs}()} first in order to obtain such data.}
  \item{data}{as \eqn{u}, an \eqn{n\times d}{n x d}-matrix of data.  For
    \code{method} being \code{"mpl"}, \code{"ml"} or \code{"itau.mpl"},
    this has to be data in \eqn{[0,1]^d}.  For \code{method} being
    \code{"itau"} or \code{"irho"}, it can either be data in \eqn{[0,1]^d}
    or in the whole \eqn{d}-dimensional space.}
  \item{copula}{a \code{"\linkS4class{copula}"} object.}
  \item{error}{(for \code{loglikCopula()}:) a \code{\link{character}}
    string specifying how errors in the underlying \code{\link{dCopula}()}
    calls should be handled:
    \describe{
      \item{\code{"-Inf"}:}{the value of the log likelihood should silently be set
	to \code{-Inf}.}

      \item{\code{"warn-Inf"}:}{signal a \code{\link{warning}} about the error and
	set the value to \code{-Inf}.}
      \item{\code{"let-it-be"}:}{the error is signalled and hence the likelihood
	computation fails.}
    }
  }
  \item{method}{a \code{\link{character}} string specifying the copula
    parameter estimator used.  This can be one of:
    \describe{
      \item{"mpl"}{Maximum pseudo-likelihood estimator (based on
	\dQuote{pseudo-observations} in \eqn{[0,1]^d}, typical obtained
	via \code{\link{pobs}()}).}
      \item{"ml"}{As \code{"mpl"} just with a different variance
	estimator.  For this to be correct (thus giving the true MLE),
	\code{data} are assumed to be observations from the true
	underlying copula whose parameter is to be estimated.}
      \item{"itau"}{Inversion of Kendall's tau estimator. \code{data}
	can be either in \eqn{[0,1]^d} (true or pseudo-observations of
	the underlying copula to be estimated) or in the \eqn{d}-dimensional space.}
      \item{"irho"}{As \code{"itau"} just with Spearman's rho instead of
	Kendall's tau.}
      \item{"itau.mpl"}{This is the estimator of \eqn{t} copula
	parameters suggested by Mashal and Zeevi (2002) based on the
	idea of inverting Kendall's tau for estimating the correlation
	matrix as introduced in a RiskLab report in 2001 later published
	as Embrechts et al. (2003); see also
	Demarta and McNeil (2005). The
	given \code{data} has to be in \eqn{[0,1]^d} (either true or
	pseudo-observations of the underlying copula to be estimated).
	Note that this method requires \code{dispstr = "un"}.}
    }
  }
  \item{posDef}{a \code{\link{logical}} indicating whether a proper
    correlation matrix is computed.}
  \item{start}{a \code{\link{vector}} of starting values for the
    parameter optimization via \code{\link{optim}()}.}
  \item{lower, upper}{Lower or upper parameter bounds for the
    optimization methods \code{"Brent"} or \code{"L-BFGS-B"}.}
  \item{optim.control}{a \code{\link{list}} of control parameters
    passed to \code{\link{optim}(*, control=optim.control)}.}
  \item{optim.method}{a character string specify the optimization
    method \emph{or} a \code{\link{function}} which when called with
    arguments \code{(copula, method, dim)} will return such a character string,
    see \code{\link{optim}()}'s \code{method}; only used when
    \code{method = "mpl"} or \code{"ml"}.

    The default has been changed (for \pkg{copula} 0.999-16, in
    Aug. 2016) from \code{"BFGS"} to the result of
    \code{optimMeth(copula, method, dim)} which is often
    \code{"L-BFGS-B"}.}
  \item{dim}{integer, the data and copula dimension, \eqn{d \ge 2}{d >= 2}.}

  \item{estimate.variance}{a \code{\link{logical}} indicating whether
    the estimator's asymptotic variance is computed (if available for
    the given \code{copula}; the default \code{NA} computes it for the
    \code{method}s \code{"itau"} and \code{"irho"}, cannot (yet) compute
    it for \code{"itau.mpl"} and only computes it for \code{"mpl"} or \code{"ml"}
    if the optimization converged).}
  \item{hideWarnings}{a \code{\link{logical}}, which, if
    \code{\link{TRUE}}, suppresses warnings from the involved
    likelihood maximization (typically when the likelihood is
    evaluated at invalid parameter values).}
  \item{\dots}{additional arguments passed to \code{method} specific auxiliary
    functions, e.g., \code{traceOpt = TRUE} (or \code{traceOpt = 10}) for tracing
    \code{\link{optimize}} (every 10-th function evaluation) for method
    \code{"itau.mpl"}, and for \dQuote{"manual"} tracing with method
    \code{"ml"} or \code{"mpl"} (notably for \code{optim.method="Brent"}).}
}
\value{
  \code{loglikCopula()} returns the copula log-likelihood evaluated at the
  parameter (vector) \code{param} given the data \code{u}.

  The return value of \code{fitCopula()} is an object of class
  \code{"\linkS4class{fitCopula}"} (inheriting from hidden class
  \code{"fittedMV"}), containing (among others!) the slots
  \describe{
    \item{estimate}{The parameter estimates.}
    \item{var.est}{The large-sample (i.e., asymptotic) variance estimate of the parameter
      estimator; (filled with) \code{NA} if \code{estimate.variance=FALSE}.
    }
    \item{copula}{The fitted copula object.}
  }

  The \code{\link{summary}()} method for \code{"fitCopula"} objects
  returns an S3 \dQuote{class} \code{"summary.fitCopula"}, which is simply a list
  with components \code{method}, \code{loglik} and \code{convergence},
  all three from the corresponding slots of the
  \code{"\linkS4class{fitCopula}"} objects, and \code{coefficients}
  (a matrix of estimated coefficients, standard errors, t values and p-values).
}
\details{
  The only difference between \code{"mpl"} and \code{"ml"} is in the
  variance-covariance estimate, \emph{not} in the parameter
  (\eqn{\theta}) estimates.

  If method \code{"mpl"} in \code{fitCopula()} is used and if
  \code{start} is not assigned a value, estimates obtained from method
  \code{"itau"} are used as initial values in the optimization. Standard
  errors are computed as explained in Genest, Ghoudi and Rivest (1995);
  see also Kojadinovic and Yan (2010, Section 3).  Their estimation
  requires the computation of certain partial derivatives of the (log)
  density.  These have been implemented for six copula families thus far:
  the Clayton, Gumbel-Hougaard, Frank, Plackett, normal and \eqn{t} copula
  families.  For other families, numerical differentiation based on
  \code{\link[numDeriv]{grad}()} from package \CRANpkg{numDeriv}
  is used (and a warning message is displayed).

  In the multiparameter elliptical case and when the estimation is based
  on Kendall's tau or Spearman's rho, the estimated correlation matrix
  may not always be positive-definite. In that case,
  \code{\link[Matrix]{nearPD}(*, corr=TRUE)} (from \CRANpkg{Matrix}) is
  applied to get a proper correlation matrix.

  For normal and \eqn{t} copulas, \code{fitCopula(, method = "mpl")} and
  \code{fitCopula(, method = "ml")} maximize the log-likelihood based on
  \pkg{mvtnorm}'s \code{dmvnorm()} and \code{dmvt()}, respectively.
  The latter two functions set the respective densities to zero if the
  correlation matrices of the corresponding distributions are not
  positive definite. As such, the estimated correlation matrices
  will be positive definite.

  If methods \code{"itau"} or \code{"irho"} are used in
  \code{fitCopula()}, an estimate of the asymptotic variance (if
  available for the copula under consideration) will be correctly
  computed only if the argument \code{data} consists of pseudo-observations
  (see \code{\link{pobs}()}).

  Consider the \eqn{t} copula with \code{df.fixed=FALSE} (see
  \code{\link{ellipCopula}()}).  In this case, the methods \code{"itau"}
  and \code{"irho"} cannot be used in \code{fitCopula()} as they cannot
  estimate the degrees of freedom parameter \code{df}.  For the methods
  \code{"mpl"} and \code{"itau.mpl"} the asymptotic variance cannot be
  (fully) estimated (yet).  For the methods \code{"ml"} and \code{"mpl"},
  when \code{start} is not specified, the starting value for \code{df}
  is set to \code{copula@df}, typically 4.

  To implement the \emph{Inference Functions for Margins} (IFM) method
  (see, e.g., Joe 2005), set \code{method="ml"} and note that
  \code{data} need to be parametric pseudo-observations obtained from \emph{fitted}
  parametric marginal distribution functions.  The returned large-sample
  variance will then underestimate the true variance (as the procedure
  cannot take into account the (unknown) estimation error for the margins).

  The fitting procedures based on \code{\link{optim}()} generate
  warnings because invalid parameter values are tried during the
  optimization process.  When the number of parameters is one and the
  parameter space is bounded, using \code{optim.method="Brent"} is
  likely to give less warnings.  Furthermore, from experience,
  \code{optim.method="Nelder-Mead"} is sometimes a more robust
  alternative to \code{optim.method="BFGS"} or \code{"L-BFGS-B"}.

  There are methods for \code{\link{vcov}()}, \code{\link{coef}()},
  \code{\link{logLik}()}, and \code{\link{nobs}()}.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Genest, C. (1987). Frank's family of bivariate distributions.
  \emph{Biometrika} \bold{74}, 549--555.

  Genest, C. and Rivest, L.-P. (1993). Statistical inference procedures
  for bivariate Archimedean copulas.
  \emph{Journal of the American Statistical Association} \bold{88}, 1034--1043.

  Rousseeuw, P. and Molenberghs, G. (1993). Transformation of nonpositive
  semidefinite correlation matrices.
  \emph{Communications in Statistics: Theory and Methods} \bold{22}, 965--984.

  Genest, C., Ghoudi, K., and Rivest, L.-P. (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika} \bold{82}, 543--552.

  Joe, H. (2005). Asymptotic efficiency of the two-stage estimation
  method for copula-based models.
  \emph{Journal of Multivariate Analysis} \bold{94}, 401--419.

  Mashal, R. and Zeevi, A. (2002). Beyond Correlation: Extreme Co-movements Between
  Financial Assets.
  \url{https://www0.gsb.columbia.edu/faculty/azeevi/PAPERS/BeyondCorrelation.pdf} (2016-04-05)

  Demarta, S. and McNeil, A. J. (2005). The t copula and related copulas.
  \emph{International Statistical Review} \bold{73}, 111--129.

  Genest, C. and Favre, A.-C. (2007). Everything you always wanted to know about
  copula modeling but were afraid to ask.
  \emph{Journal of Hydrologic Engineering} \bold{12}, 347--368.

  Kojadinovic, I. and Yan, J. (2010). Comparison of three semiparametric methods
  for estimating dependence parameters in copula models.
  \emph{Insurance: Mathematics and Economics} \bold{47}, 52--63.
}
\seealso{
  \code{\link{Copula}},
  \code{\link{fitMvdc}} for fitting multivariate distributions
  \emph{including} the margins,
  \code{\link{gofCopula}} for goodness-of-fit tests.

  For maximum likelihood of (nested) Archimedean copulas, see
  \code{\link{emle}}, etc.
}
\examples{
(Xtras <- copula:::doExtras()) # determine whether examples will be extra (long)
n <- if(Xtras) 200 else 64 # sample size

## A Gumbel copula
set.seed(7) # for reproducibility
gumbel.cop <- gumbelCopula(3, dim=2)
x <- rCopula(n, gumbel.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(gumbelCopula(), u, method="itau")
fit.tau
confint(fit.tau) # work fine !
confint(fit.tau, level = 0.98)
summary(fit.tau) # a bit more, notably "Std. Error"s
coef(fit.tau)# named vector
coef(fit.tau, SE = TRUE)# matrix
\dontshow{stopifnot(is.matrix(cm <- coef(fit.tau, SE = TRUE)), cm[,2] > 0)}
## Inverting Spearman's rho
fit.rho <- fitCopula(gumbelCopula(), u, method="irho")
summary(fit.rho)
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(gumbelCopula(), u, method="mpl")
fit.mpl
## Maximum likelihood -- use 'x', not 'u' ! --
fit.ml <- fitCopula(gumbelCopula(), x, method="ml")
summary(fit.ml) # now prints a bit more than simple 'fit.ml'
## ... and what's the log likelihood (in two different ways):
(ll. <- logLik(fit.ml))
stopifnot(all.equal(as.numeric(ll.),
            loglikCopula(coef(fit.ml), u=x, copula=gumbel.cop)))

## A Gauss/normal copula

## With multiple/*un*constrained parameters
set.seed(6) # for reproducibility
normal.cop <- normalCopula(c(0.6, 0.36, 0.6), dim=3, dispstr="un")
x <- rCopula(n, normal.cop) # "true" observations (simulated)
u <- pobs(x)                # pseudo-observations
## Inverting Kendall's tau
fit.tau <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="irho")
fit.rho
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="un"), u, method="mpl")
summary(fit.mpl)
coef(fit.mpl) # named vector
coef(fit.mpl, SE = TRUE) # the matrix, with SE
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="un"), x, method="ml", traceOpt=TRUE)
summary(fit.ml)
confint(fit.ml)
confint(fit.ml, level = 0.999) # clearly non-0

## Fix some of the parameters
param <- c(.6, .3, NA_real_)
attr(param, "fixed") <- c(TRUE, FALSE, FALSE)
ncp <- normalCopula(param = param, dim = 3, dispstr = "un")
fixedParam(ncp) <- c(TRUE, TRUE, FALSE)
## 'traceOpt = 5': showing every 5-th log likelihood evaluation:
summary(Fxf.mpl <- fitCopula(ncp, u, method = "mpl", traceOpt = 5))
Fxf.mpl@copula # reminding of the fixed param. values

## With dispstr = "toep" :
normal.cop.toep <- normalCopula(c(0, 0), dim=3, dispstr="toep")
## Inverting Kendall's tau
fit.tau <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="itau")
fit.tau
## Inverting Spearman's rho
fit.rho <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="irho")
summary(fit.rho)
## Maximum pseudo-likelihood
fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="toep"), u, method="mpl")
fit.mpl
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="toep"), x, method="ml")
summary(fit.ml)

## With dispstr = "ar1"
normal.cop.ar1 <- normalCopula(c(0), dim=3, dispstr="ar1")
## Inverting Kendall's tau
summary(fit.tau <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="itau"))
## Inverting Spearman's rho
summary(fit.rho <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="irho"))
## Maximum pseudo-likelihood
summary(fit.mpl <- fitCopula(normalCopula(dim=3, dispstr="ar1"), u, method="mpl"))
## Maximum likelihood (use 'x', not 'u' !)
fit.ml <- fitCopula(normalCopula(dim=3, dispstr="ar1"), x, method="ml")
summary(fit.ml)

## A t copula with variable df (df.fixed=FALSE)
(tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5))
set.seed(101)
x <- rCopula(n, tCop) # "true" observations (simulated)
## Maximum likelihood (start = (rho[1:3], df))
summary(tc.ml <- fitCopula(tCopula(dim=3, dispstr="un"), x, method="ml",
                           start = c(0,0,0, 10)))
## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
u <- pobs(x)          # pseudo-observations
tc.mpl <- fitCopula(tCopula(dim=3, dispstr="un"),
                     u, method="mpl", estimate.variance=FALSE,
                     start= c(0,0,0, 10))
summary(tc.mpl)
}% extra checks etc ---> ../tests/fitting-ex.R
\keyword{models}
\keyword{multivariate}
