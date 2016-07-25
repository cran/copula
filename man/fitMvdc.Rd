\name{fitMvdc}
\title{Estimation of Multivariate Models Defined via Copulas}
\alias{fitMvdc}
\alias{coef.fittedMV}
\alias{logLik.fittedMV}
\alias{vcov.fittedMV}
\alias{loglikMvdc}
\description{
  Fitting copula-based multivariate distributions
  (\code{"\linkS4class{mvdc}"}) to multivariate data,
  estimating both the marginal and the copula parameters.

  If you assume {non parametric} margins, in other words, take the
  empirical distributions for all margins, you can use
  \code{\link{fitCopula}(*, \link{pobs}(x))} instead.
}
\usage{
loglikMvdc(param, x, mvdc)
fitMvdc(data, mvdc, start, optim.control = list(), method = "BFGS",
        lower = -Inf, upper = Inf,
        estimate.variance = fit$convergence == 0, hideWarnings = TRUE)

\S3method{coef}{fittedMV}(object, ...)
\S3method{logLik}{fittedMV}(object, ...)
\S3method{vcov}{fittedMV}(object, ...)
}
\arguments{
  \item{param}{a vector of parameter values.  When specifying parameters for
    \code{mvdc} objects, the parameters must be ordered with the
    marginals first and the copula parameters last.  When the
    \code{mvdc} object has \code{marginsIdentical = TRUE}, only the
    parameters of one marginal must be specified.}
  \item{x}{a data matrix.}
  \item{mvdc}{a \code{"\linkS4class{mvdc}"} object. }
  \item{hideWarnings}{logical indicating if warning messages from
    likelihood maximization, e.g., from evaluating at invalid parameter
    values, should be suppressed (via \code{\link{suppressWarnings}}).}
  \item{data}{a data matrix.}
  \item{start}{a vector of starting value for \code{"param"}.  See
    \code{"param"} above for ordering of this vector.}
  \item{optim.control}{a list of controls to be passed to \code{\link{optim}}.}
  \item{method}{the method for \code{\link{optim}}.}
  \item{lower, upper}{bounds on each parameter, passed to
    \code{\link{optim}}, typically \dQuote{box constraints} for
    \code{method = "L-BFGS-B"}.}
  \item{estimate.variance}{logical; if true (as by default, if the
    optimization converges), the asymptotic variance is estimated.}

  \item{object}{an \R object of class \code{"\linkS4class{fitMvdc}"}.}
  \item{\dots}{potentially further arguments to methods.}
}
\value{
  The return value \code{loglikMvdc()} is the log likelihood evaluated
  for the given value of \code{param}.

  The return value of \code{fitMvdc()} is an object of class
  \code{"\linkS4class{fitMvdc}"} % -> ./fitCopula-class.Rd
  (see there), containing slots (among others!):
  \item{estimate}{the estimate of the parameters.}
  \item{var.est}{large-sample (i.e., asymptotic) variance estimate of the parameter
    estimator (filled with \code{NA} if \code{estimate.variance = FALSE}).}
  \item{mvdc}{the \emph{fitted} multivariate distribution, see
    \code{\linkS4class{mvdc}}.}

%% FIXME -- formatting -- because of the special way \item is handled
%% inside \value{...}, there is no empty line put here, even though we
%% want it.

  The \code{\link{summary}()} method for \code{"fitMvdc"} objects
  returns a S3 \dQuote{class} \code{"summary.fitMvdc"}, simply a list
  with components \code{method}, \code{loglik}, and \code{convergence},
  all three from corresponding slots of the
  \code{"\linkS4class{fitMvdc}"} objects, and
  \item{coefficients}{a matrix of estimated coefficients, standard
  errors, t values and p-values.}

}
\note{
  User-defined marginal distributions can be used as long as the
  \code{"{dpq}"} functions are defined.
  See \code{vignette("AR_Clayton", package="copula")}.

  When covariates are available for marginal distributions or for the copula,
  one can construct the loglikelihood function and feed it to \code{"optim"}
  to estimate all the parameters.

  Finally, note that some of the fitting functions generate error
  messages because invalid parameter values are tried during the
  optimization process (see \code{\link{optim}}).  This should be rarer
  since 2013, notably for likelihood based methods (as the likelihood
  is now rather set to \code{-Inf} than giving an error).

  Previously, \code{loglikMvdc()} had an argument \code{hideWarnings};
  nowadays, do use \code{\link{suppressWarnings}(..)} if you are sure
  you do not want to see them.
}
%\author{Ivan Kojadinovic and Martin Maechler}
\seealso{
  \code{\link{mvdc}} and \code{\linkS4class{mvdc}};
  further, \code{\link{Copula}}, \code{\link{fitCopula}},
  \code{\link{gofCopula}}.

  For fitting univariate marginals, \code{\link[MASS]{fitdistr}()}.
}
\examples{
gumbel.cop <- gumbelCopula(3, dim=2)
gMvd2 <- mvdc(gumbel.cop, c("exp","exp"),
              list(list(rate=2), list(rate=4)))
set.seed(11)
x <- rMvdc(10000, gMvd2)
## with identical margins:
gMvd.I <- mvdc(gumbel.cop, "exp", param = list(rate=3), marginsIdentical=TRUE)
if(copula:::doExtras()) { ## these are typically not run with CRAN checking:
 ## Takes about 25 sec. [2012-07]:
 fit2 <- fitMvdc(x, gMvd2, start = c(1,1, 2),
                 hideWarnings=FALSE) ## <- show warnings here
 print(fit2)
 ## The (estimated, asymptotic) var-cov matrix:
 print( vcov(fit2) )

 fitI <- fitMvdc(x, gMvd.I, start = c(3, 2),
                 optim.control=list(trace= TRUE, REPORT= 2))
 print(coef(summary(fitI)))
 print(fitI)

  ## a wrong starting value can already be *the* problem:
  f2 <- try(fitMvdc(x, gMvd.I, start = c(1, 1),
	            optim.control=list(trace= TRUE, REPORT= 2)))
  ##--> Error in optim( ... ) : non-finite finite-difference value [2]

  ##==> "Solution" :  Using a more robust (but slower) optim() method:
  fI.2 <- fitMvdc(x, gMvd.I, start = c(1, 1), method = "Nelder",
		  optim.control=list(trace= TRUE))

 ## The (estimated, asymptotic) var-cov matrix:
 print( vcov(fit2) )
 str(sfI <- summary(fitI))
 stopifnot(is.matrix(coef(sfI)))
\dontshow{
 stopifnot(all.equal(coef(sfI),
		     coef(summary(fI.2)),    tolerance = 2e-4),
	   all.equal(coef(fitI), coef(fI.2), tolerance = 1e-4))
}%dont...
}% if(..doExtras())
}
\keyword{models}
\keyword{multivariate}
