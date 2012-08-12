\name{fitCopula}
\alias{fitCopula}
\alias{loglikCopula}
\title{Estimation of the parameters in copula models}
\description{
  Fits a copula model to multivariate data belonging to the unit
  hypercube. The data can be pseudo-observations constructed from
  empirical or parametric marginal c.d.f.s, or true observations
  from the copula.
}
\usage{
loglikCopula(param, x, copula, hideWarnings = FALSE)
fitCopula(copula, data, method = c("mpl", "ml", "itau", "irho"),
          start = NULL, lower = NULL, upper = NULL,
          optim.method = "BFGS", optim.control = list(maxit=1000),
          estimate.variance = TRUE, hideWarnings = TRUE)
}
\arguments{
  \item{param}{ a vector of parameter values. }
  \item{x, data}{numeric data matrix.}
  \item{copula}{ a \code{"\linkS4class{copula}"} object. }
  \item{hideWarnings}{logical; if TRUE, warning messages from
    evaluating the loglikelihood at invalid parameter values are suppressed.}
  \item{method}{a \code{\link{character}} string specifying the method;
    can be either
    \code{"ml"} (maximum likelihood),
    \code{"mpl"} (maximum pseudo-likelihood),
    \code{"itau"} (inversion of Kendall's tau), and
    \code{"irho"} (inversion of Spearman's rho).  The last three methods
    assume that the \code{data} are pseudo-observations (scaled
    ranks), while the first method assumes that the \code{data} are
    observations from the unknown copula.  The default is \code{"mpl"}.}
  \item{start}{a vector of starting values for \code{param}.}
  \item{lower, upper}{bounds on the variables for the \code{"Brent"}
    or \code{"L-BFGS-B"} method.}
  \item{optim.control}{a \code{\link{list}} of control parameters to be
    passed to \code{\link{optim}(*, control=optim.control)}.}
  \item{optim.method}{the \code{method} for \code{\link{optim}}.}
  \item{estimate.variance}{logical; if TRUE, the asymptotic variance is
    estimated.}
}
\value{
  \code{loglikCopula()} returns the log likelihood evaluated at the
  given value of \code{"param"}.

  The return value of \code{fitCopula()} is an object of class
  \code{"\linkS4class{fitCopula}"} containing slots (among others):
  \item{estimate}{ the estimate of the parameters. }
  \item{var.est}{ estimate of the asymptotic variance estimate of the
    parameter estimator. }
  \item{method}{ the estimation method. }
  \item{loglik}{ loglikelihood at \code{"estimate"}. }
  \item{copula}{ the fitted copula. }
}
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

  Demarta, S. and McNeil, A. J. (2005). The t copula and related copulas.
  \emph{International Statistical Review} \bold{73}, 111--129.

  Genest, C. and Favre, A.-C. (2007). Everything you always wanted to know about
  copula modeling but were afraid to ask.
  \emph{Journal of Hydrologic Engineering} \bold{12}, 347--368.

  Kojadinovic, I. and Yan, J. (2010). Comparison of three semiparametric methods
  for estimating dependence parameters in copula models.
  \emph{Insurance: Mathematics and Economics} \bold{47}, 52--63.
}
\note{
  In the multiparameter elliptical case and when the estimation is based
  on Kendall's tau or Spearman's rho, the
  estimated correlation matrix may not always be positive-definite.  If
  it is not, the correction proposed by Rousseeuw and Molenberghs (1993)
  is applied and a warning message given.

  If method \code{"mpl"} in \code{fitCopula()} is used and if
  \code{start} is not assigned a value, estimates obtained
  from method \code{"itau"} are used as initial values in
  the optimization.

  For the t copula with \code{df.fixed=FALSE} (see \code{\link{ellipCopula}}),
  the methods \code{"itau"} and \code{"irho"} cannot be used in
  \code{fitCopula()}.  For the methods \code{"ml"} and \code{"mpl"},
  when \code{start} is not specified, the starting value for \code{df} is
  set to \code{copula@df}, typically 4.  Also, the asymptotic variance
  cannot (yet) be estimated for method \code{"mpl"}.

  To implement the \dQuote{\emph{inference functions for margins}}
  method (see e.g., Joe 2005), the \code{data} need to be pseudo-observations
  obtained from fitted parametric marginal c.d.f.s and \code{method}
  needs to be set to \code{"ml"}.  The returned large-sample variance
  will then underestimate the true variance.

  Finally, note that the fitting functions generate error messages
  because invalid parameter values are tried during the optimization
  process (see \code{\link{optim}}).  When the number of parameters is
  one and the parameter space is bounded, using
  \code{optim.method="Brent"} is likely to give less warnings.
  Furthermore, from experience, \code{optim.method="Nelder-Mead"} is
  sometimes a more robust alternative to \code{optim.method="BFGS"}.
}
\seealso{\code{\link{Copula}},
  \code{\link{mvdc}} for fitting multivariate distributions
  \emph{including} the margins aka \dQuote{meta copula}s;
  \code{\link{gofCopula}}.
}
\examples{
gumbel.cop <- gumbelCopula(3, dim=2)

n <- 200
x <- rCopula(n, gumbel.cop)## "true" observations
u <- pobs(x)               ## pseudo-observations
## inverting Kendall's tau
fit.tau <- fitCopula(gumbel.cop, u, method="itau")
fit.tau
## inverting Spearman's rho
fit.rho <- fitCopula(gumbel.cop, u, method="irho")
fit.rho
## maximum pseudo-likelihood
fit.mpl <- fitCopula(gumbel.cop, u, method="mpl")
fit.mpl
## maximum likelihood
fit.ml <- fitCopula(gumbel.cop, x, method="ml")
fit.ml

## a multiparameter example
normal.cop <- normalCopula(c(0.6,0.36, 0.6),dim=3,dispstr="un")
x <- rCopula(n, normal.cop) ## "true" observations
u <- pobs(x)                ## pseudo-observations
## inverting Kendall's tau
fit.tau <- fitCopula(normal.cop, u, method="itau")
fit.tau
## inverting Spearman's rho
fit.rho <- fitCopula(normal.cop, u, method="irho")
fit.rho
## maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop, u, method="mpl")
fit.mpl
## maximum likelihood
fit.ml <- fitCopula(normal.cop, x, method="ml")
fit.ml
## with dispstr="toep"
normal.cop.toep <- normalCopula(c(0, 0), dim=3, dispstr="toep")
## inverting Kendall's tau
fit.tau <- fitCopula(normal.cop.toep, u, method="itau")
fit.tau
## inverting Spearman's rho
fit.rho <- fitCopula(normal.cop.toep, u, method="irho")
fit.rho
## maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop.toep, u, method="mpl")
fit.mpl
## maximum likelihood
fit.ml <- fitCopula(normal.cop.toep, x, method="ml")
fit.ml
## with dispstr="ar1"
normal.cop.ar1 <- normalCopula(c(0), dim=3, dispstr="ar1")
## inverting Kendall's tau
fit.tau <- fitCopula(normal.cop.ar1, u, method="itau")
fit.tau
## inverting Spearman's rho
fit.rho <- fitCopula(normal.cop.ar1, u, method="irho")
fit.rho
## maximum pseudo-likelihood
fit.mpl <- fitCopula(normal.cop.ar1, u, method="mpl")
fit.mpl
## maximum likelihood
fit.ml <- fitCopula(normal.cop.ar1, x, method="ml")
fit.ml

## a t copula with variable df  (df.fixed=FALSE):
(tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5))
set.seed(101)
x <- rCopula(n, tCop) ## "true" observations
u <- pobs(x)          ## pseudo-observations
## maximum likelihood;   start := (rho[1:3], df)
(tc.ml  <- fitCopula(tCop, x, method="ml", start=c(0,0,0, 10)))
(tc.ml. <- fitCopula(tCop, x, method="ml")) # without 'start'
## maximum pseudo-likelihood; the asymptotic variance cannot be estimated
(tc.mpl <- fitCopula(tCop, u, method="mpl", estimate.variance=FALSE,
                     start= c(0,0,0,10)))
## without start:
(tc.mp. <- fitCopula(tCop, u, method="mpl", estimate.variance=FALSE))
stopifnot(all.equal(tc.ml , tc.ml., tol= .005),
          all.equal(tc.mpl, tc.mp., tol= .005))

## same t copula but with df.fixed=TRUE	 (--> use same data!)
(tC.f <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5, df.fixed=TRUE))
## maximum likelihood;  start := rho[1:3]		    -------------
(tcF.ml  <- fitCopula(tC.f, x, method="ml", start=c(0,0,0)))
(tcF.ml. <- fitCopula(tC.f, x, method="ml"))# without 'start'
stopifnot(all.equal(tcF.ml,tcF.ml., tol= 4e-4))

## maximum pseudo-likelihood; the asymptotic variance cannot be estimated
(tcF.mpl <- fitCopula(tC.f, u, method="mpl", estimate.variance=FALSE,
                      start=c(0,0,0)))
(tcF.mp. <- fitCopula(tC.f, u, method="mpl", estimate.variance=FALSE))
stopifnot(all.equal(tcF.mpl,tcF.mp., tol= 1e-5))
}
\keyword{ models }
\keyword{ multivariate }
