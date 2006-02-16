\name{fitCopula}
\alias{fitCopula}
\alias{fitMvdc}
\alias{loglikCopula}
\alias{loglikMvdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum Likelihood Estimation of Copula Models }
\description{
  Fit a copula model to multivariate data.
}
\usage{
loglikCopula(param, x, copula)
loglikMvdc(param, x, mvdc)
fitCopula(data, copula, start, optim.control = list(NULL), method = "BFGS")
fitMvdc(data, mvdc, start, optim.control = list(NULL), method = "BFGS")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{param}{ a vector of parameter values }
  \item{x}{ a data matrix }
  \item{copula}{ a 'copula' object }
  \item{mvdc}{ a 'mvdc' object }

  \item{data}{ a data matrix }
  \item{start}{ a vector of starting value for param }
  \item{optim.control}{ a list of control to be passed to optim }
  \item{method}{ the method for optim }
}
\value{
  The return values of 'loglikCopula' and 'loglikMvdc' are the
  loglikelihood evaluated at the given value of 'param'.
  
  The return values of 'fitCopula' and 'fitMvdc' are a list of:
  \item{est}{ the estimate of the parameters }
  \item{var.est}{ variance matrix of the estimate }
  \item{loglik}{ loglikelihood at est }
  \item{fit}{ the result of optim }
}
\references{
  Yan (2006) Multivariate Modeling with Copulas and Engineering
  Applications. In Handbook of Engineering Statistics, Ed. Pham, Springer. 
  }
\author{ Jun Yan <jyan@stat.uiowa.edu>}
\note{
  When covariates are available for marginal distributions or copula,
  one can construct loglikelihood function and feed it to optim to
  estimate all the parameters.
}
\seealso{ \code{\link{Copula}}, \code{\link{mvdc}} }
\examples{
gmb <- gumbelCopula(3, dim=2)
myMvd <- mvdc(gmb, c("exp","exp"), list(list(rate=2),list(rate=4)))
x <- rmvdc(myMvd, 1000)
fit <- fitMvdc(x, myMvd, c(1,1,2))
}
\keyword{ models}% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
