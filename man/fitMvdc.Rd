\name{fitMvdc}
\alias{fitMvdc}
\alias{loglikMvdc}
\title{ Estimation of multivariate models defined via copulas }
\description{
  Fits a copula-based multivariate distribution to multivariate data.
}
\usage{
loglikMvdc(param, x, mvdc, suppressMessages=FALSE)
fitMvdc(data, mvdc, start, optim.control = list(), method = "BFGS")
}
\arguments{
  \item{param}{ a vector of parameter values.  When specifying parameters for
    \code{mvdc} objects the parameters should be ordered with the marginals first
    and the copula parameters last.  When the \code{mvdc} object has
    \code{marginsIdentical = TRUE}, only the parameters of one marginal
    should be specified.}
  \item{x}{ a data matrix. }
  \item{mvdc}{ a \code{"\linkS4class{mvdc}"} object. }
  \item{suppressMessages}{ logical, if TRUE, warnings messages from
    evaluating loglikelihood at invalid parameter values are
    suppressed. }
  \item{data}{ a data matrix. }
  \item{start}{ a vector of starting value for \code{"param"}. See \code{"param"}
	above for ordering of this vector.}
  \item{optim.control}{ a list of controls to be passed to \code{\link{optim}}. }
  \item{method}{ the method for \code{\link{optim}}. }
}
\value{
  The return value \code{loglikMvdc()} is the log likelihood evaluated
  for the given value of \code{param}.

  The return value of \code{fitMvdc()} is an object of class
  \code{"\linkS4class{fitMvdc}"}, % -> ./fitCopula-class.Rd
  containing slots:
  \item{estimate}{ the estimate of the parameters. }
  \item{var.est}{ large-sample variance estimate of the parameter estimator. }
  \item{loglik}{ loglikelihood at \code{"est"}. }
  \item{copula}{ the fitted copula. }
}

\note{
  User-defined marginal distribution can be used as long as the
  \code{"{dpq}"} functions are defined.  See \code{demo(QARClayton)}
  prepared by Roger Koenker \email{rkoenker@uiuc.edu}.

  When covariates are available for marginal distributions or for the copula,
  one can construct the loglikelihood function and feed it to \code{"optim"}
  to estimate all the parameters.

  Finally, note that the fitting functions generate error messages
  because invalid parameter values are tried during the optimization
  process (see \code{\link{optim}}).
}
\seealso{ \code{\link{Copula}}, \code{\link{fitCopula}}, \code{\link{gofCopula}}. }

\examples{
gumbel.cop <- gumbelCopula(3, dim=2)
myMvd <- mvdc(gumbel.cop, c("exp","exp"), list(list(rate=2),list(rate=4)))
x <- rmvdc(myMvd, 1000)
fit <- fitMvdc(x, myMvd, c(1,1,2))
fit

## Roger Koenker prepared a demo illustrating MLE for a Clayton AR(1)
## copula model with identical, user-defined Student marginals
if(FALSE) ## not yet ok %% FIXME !
demo("QARClayton")
}
\keyword{ models }
\keyword{ multivariate }
