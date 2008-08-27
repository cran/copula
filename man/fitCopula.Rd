\name{fitCopula}
\alias{fitCopula}
\alias{fitMvdc}
\alias{loglikCopula}
\alias{loglikMvdc}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Maximum likelihood estimation of copula models }
\description{
  Fits a copula model to multivariate data.
}
\usage{
loglikCopula(param, x, copula)
loglikMvdc(param, x, mvdc)
fitCopula(data, copula, start, lower=NULL, upper=NULL,
          optim.control = list(NULL), method = "BFGS")
fitMvdc(data, mvdc, start, optim.control = list(NULL), method = "BFGS")
}
\arguments{
  \item{param}{ a vector of parameter values. }
  \item{x}{ a data matrix. }
  \item{copula}{ a \code{"copula"} object. }
  \item{mvdc}{ a \code{"mvdc"} object. }

  \item{data}{ a data matrix. }
  \item{start}{ a vector of starting value for \code{param}. }
  \item{lower, upper}{ bounds on the variables for the \code{"L-BFGS-B"} method.}
  \item{optim.control}{ a list of controls to be passed to \code{"optim"}. }
  \item{method}{ the method for optim. }
}
\value{
  The return value of \code{"loglikCopula"} and \code{"loglikMvdc"} is
  the loglikelihood evaluated for the given value of \code{"param"}.
    
  The return values of \code{"fitCopula"} and \code{"fitMvdc"} are
  objects of class \code{"fitCopula"} and \code{"fitMvdc"}, respectively,
  containing slots:
  \item{est}{ the estimate of the parameters. }
  \item{var.est}{ variance matrix estimate of
      the parameter estimator. See note below. }
  \item{loglik}{ loglikelihood at \code{"est"}. }
  \item{fit}{ the result of \code{"optim"}. }
}
\references{
  C. Genest, K. Ghoudi and L.-P. Rivest (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika}, 82, 543-552. 

  J. Yan (2006), Multivariate Modeling with Copulas and Engineering
  Applications. In \emph{Handbook of Engineering Statistics}, Ed. Pham,
  Springer.
}
\note{
  The variance matrix of the estimate, \code{"var.est"}, in the returned
  object is only valid for full likelihood maximization. When
  pseudo-observations are used, they under-estimate the variation. See
  Genest, Ghoudi, and Rivest (1995).
  
  When covariates are available for marginal distributions or for the copula,
  one can construct the loglikelihood function and feed it to \code{"optim"}
  to estimate all the parameters.
}

\seealso{ \code{\link{Copula}}, \code{\link{mvdc}}, \code{\link{gofCopula}}. }

\examples{
gumbel.cop <- gumbelCopula(3, dim=2)
myMvd <- mvdc(gumbel.cop, c("exp","exp"), list(list(rate=2),list(rate=4)))
x <- rmvdc(myMvd, 1000)
fit <- fitMvdc(x, myMvd, c(1,1,2))
fit

n <- 200
u <- rcopula(gumbel.cop, n)              ## true observations
r <- apply(u, 2, rank) / (n + 1)  ## pseudo-observations
start <- calibKendallsTau(gumbel.cop, cor(u, method="kendall")[1,2])
fit.u <- fitCopula(u, gumbel.cop, start)
fit.u
fit.r <- fitCopula(r, gumbel.cop, start)
fit.r
}
\keyword{ models}% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
