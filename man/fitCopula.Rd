\name{fitCopula}
\alias{fitCopula}
\alias{loglikCopula}
\title{ Estimation of the parameters in copula models }
\description{
  Fits a copula model to multivariate data belonging to the unit
  hypercube. The data can be pseudo-observations constructed from
  empirical or parametric marginal c.d.f.s, or true observations
  from the copula. 
}
\usage{
loglikCopula(param, x, copula, hideWarnings=FALSE)
fitCopula(copula, data, method="mpl", start=NULL, lower=NULL,
upper=NULL, optim.control=list(maxit=1000), optim.method="BFGS",
estimate.variance=TRUE, hideWarnings=TRUE)
}
\arguments{
  \item{param}{ a vector of parameter values. }
  \item{x}{ a data matrix. }
  \item{data}{ a data matrix. }
  \item{copula}{ a \code{"copula"} object. }
  \item{hideWarnings}{ logical; if TRUE, warnings messages from
    evaluating loglikelihood at invalid parameter values are
    suppressed. }
  \item{method}{ can be either \code{"ml"} (maximum likelihood), \code{"mpl"} (maximum
    pseudo-likelihood), \code{"itau"} (inversion of Kendall's tau), and
    \code{"irho"} (inversion of Spearman's rho). The last three methods
    assume that the \code{"data"} are pseudo-observations (scaled
    ranks), while the first method assumes that the \code{"data"} are
    observations from the unknown copula. The default method is "mpl". }
  \item{start}{ a vector of starting value for \code{"param"}.}
  \item{lower, upper}{ bounds on the variables for the \code{"Brent"}
    or \code{"L-BFGS-B"} method.}
  \item{optim.control}{ a list of controls to be passed to \code{"optim"}. }
  \item{optim.method}{ the method for \code{"\link{optim}"}. }
  \item{estimate.variance}{ logical; if TRUE, the asymptotic variance is
    estimated. }    
}
\value{
  The return value of \code{"loglikCopula"} is
  the loglikelihood evaluated at the given value of \code{"param"}.
    
  The return values of \code{"fitCopula"} is an object of class
  \code{"fitCopula"} containing slots:
  \item{estimate}{ the estimate of the parameters. }
  \item{var.est}{ estimate of the asymptotic variance estimate of the parameter estimator. }
  \item{method}{ the estimation method. }
  \item{loglik}{ loglikelihood at \code{"estimate"}. }
  \item{copula}{ the fitted copula. }
}
\references{
  C. Genest (1987). Frank's family of bivariate
  distributions. \emph{Biometrika} 74:549--555.

  C. Genest and L.-P. Rivest (1993). Statistical inference procedures
  for bivariate Archimedean copulas.
  \emph{Journal of the American Statistical Association}
  88:1034--1043.
  
  P. Rousseeuw and G. Molenberghs (1993). Transformation of nonpositive
  semidefinite correlation matrices. \emph{Communications in Statistics: Theory
    and Methods} 22:965--984.
  
  C. Genest, K. Ghoudi and L.-P. Rivest (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika}, 82:543--552.

  H. Joe (2005). Asymptotic efficiency of the two-stage estimation
  method for copula-based models. \emph{Journal of Multivariate
    Analysis} 94:401--419.
  
  S. Demarta and A. McNeil (2005). The t copula and related copulas.
  \emph{International Statistical Review} 73:111--129.
  
  C. Genest and A.-C. Favre (2007). Everything you always wanted to know about
  copula modeling but were afraid to ask. \emph{Journal of Hydrologic Engineering}
  12:347--368.
  
  I. Kojadinovic and J. Yan (2010). Comparison of three semiparametric methods
  for estimating dependence parameters in copula models.
  \emph{Insurance: Mathematics and Economics}, 47:52--63.
}

\note{
  In the multiparameter elliptical case and when the estimation is based
  on Kendall's tau or Spearman's rho, the
  estimated correlation matrix may not always be positive-definite. If
  it is not, the correction proposed by Rousseeuw and Molenberghs (1993)
  is applied and a warning message given.

  If method \code{"mpl"} in \code{"fitCopula"} is used and if
  \code{"start"} is not assigned a value, estimates obtained
  from method \code{"itau"} are used as initial values in
  the optimization.

  For the t copula with \code{df.fixed=FALSE} (see \code{\link{ellipCopula}}),
  the methods \code{"itau"} and \code{"irho"} in \code{"fitCopula"}
  cannot be used. The methods \code{"ml"} and \code{"mpl"} can be
  employed provided the parameter \code{"start"} is assigned a
  value. See the example below. Also, the asymptotic variance cannot be
  estimated if method \code{"mpl"} is chosen.

  To implement the \emph{inference functions for margins} method (see
  e.g. Joe 2005), the \code{"data"} need to be pseudo-observations
  obtained from fitted parametric marginal c.d.f.s and \code{"method"}
  need to be set to \code{"ml"}. The returned large-sample variance
  will then underestimate the true variance.
  
  Finally, note that the fitting functions generate error messages
  because invalid parameter values are tried during the optimization
  process (see \code{\link{optim}}). When the number of parameters is
  one and the parameter space is bounded, using
  \code{optim.method="Brent"} is likely to give less warnings.
  Furthermore, from experience, \code{optim.method="Nelder-Mead"} is
  sometimes a more robust alternative to \code{optim.method="BFGS"}.
}
\seealso{ \code{\link{Copula}}, \code{\link{mvdc}}, \code{\link{gofCopula}}. }

\examples{
gumbel.cop <- gumbelCopula(3, dim=2)

n <- 200
x <- rcopula(gumbel.cop, n)       ## true observations
u <- apply(x, 2, rank) / (n + 1)  ## pseudo-observations
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
x <- rcopula(normal.cop, n)     ## true observations
u <- apply(x, 2, rank) / (n + 1)  ## pseudo-observations
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

## a t copula with df.fixed=FALSE
t.cop <- tCopula(c(0.2,0.4,0.6),dim=3,dispstr="un",df=5)
x <- rcopula(t.cop, n)     ## true observations
u <- apply(x, 2, rank) / (n + 1)  ## pseudo-observations
## maximum likelihood
fit.ml <- fitCopula(t.cop, x, method="ml", start=c(0,0,0,10))
fit.ml
## maximum pseudo-likelihood; the asymptotic variance cannot be estimated 
fit.mpl <- fitCopula(t.cop, u, method="mpl", start=c(0,0,0,10),
                     estimate.variance=FALSE)
fit.mpl
}
\keyword{ models }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
