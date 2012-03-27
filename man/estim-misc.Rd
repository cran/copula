\name{estim.misc}
\alias{ebeta}
\alias{edmle}
\alias{etau}
\title{Various Estimators for (Nested) Archimedean Copulas}
\description{
  Various Estimators for (Nested) Archimedean Copulas, namely,
  \describe{
    \item{ebeta}{Method-of-moments-like estimator based on (a
      multivariate version of) Blomqvist'sbeta.}
    \item{edmle}{Maximum likelihood estimator based on the diagonal of a
      (nested) Archimedean copula.}
    \item{etau}{Method-of-moments-like estimators based on (bivariate) Kendall's tau.}
  }
}
\usage{
ebeta(u, cop, interval = initOpt(cop@copula@name), ...)
edmle(u, cop, interval = initOpt(cop@copula@name), warn=TRUE, ...)
etau(u, cop, method = c("tau.mean", "theta.mean"), warn=TRUE, ...)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the
    sample size and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated
    (currently only Archimedean copulas are provided).}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in Hofert et
    al. (2011a).}
  \item{method}{a character string specifying the method (only
    for \code{etau}), which has to be one (or a unique abbreviation) of
    \describe{
      \item{\code{"tau.mean"}}{method-of-moments-like estimator based on
	the average of pairwise sample versions of Kendall’s tau;}
      \item{\code{"theta.mean"}}{average of the method-of-moments-like
	Kendall's tau estimators.}
    }
  }
  \item{warn}{logical indicating if warnings are printed:
    \describe{
      \item{\code{edmle()}}{for the family of \code{"Gumbel"} if the
	diagonal maximum-likelihood estimator is smaller than 1.}
      \item{\code{etau()}}{for the family of \code{"AMH"} if tau is outside
	\eqn{[0, 1/3]} and in general if at least one of the computed pairwise sample
	versions of Kendall's tau is negative.}
    }
  }
  \item{\dots}{additional arguments passed to
    \code{\link{cor}} (for \code{etau}), to \code{\link{optimize}}
    (for \code{edmle}), or to \code{\link{safeUroot}} (for \code{ebeta}).}
}
\details{
  For \code{ebeta}, the parameter is estimated with a
  method-of-moments-like procedure such that the population version of
  the multivariate Blomqvist's beta matches its sample version.

  Note that the copula diagonal is a distribution function and the
  maximum of all components of a random vector following the copula is
  distributed according to this distribution function.  For
  \code{edmle}, the parameter is estimated via maximum-likelihood
  estimation based on the diagonal.

  For \code{etau}, the \code{method="tau.mean"} means that the average
  of sample versions of Kendall's tau are computed first and then the
  parameter is determined such that the population version of Kendall's
  tau matches this average (if possible); the \code{method="theta.mean"}
  stands for first computing all pairwise Kendall's tau estimators and
  then returning the mean of these estimators.

  For more details, see Hofert et al. (2011a).

  Note that these estimators should be used with care; see the
  performance results in Hofert et al. (2011a).  In particular,
  \code{etau} should be used with the (default) method \code{"tau.mean"}
  since \code{"theta.mean"} is both slower and more prone to errors.
}
\value{
  \describe{
    \item{\code{ebeta}}{the return value of \code{\link{safeUroot}}
      (that is, typically almost the same as the value of
      \code{\link{uniroot}}) giving the Blomqvist beta estimator.}
    \item{\code{edmle}}{\code{\link{list}} as returned by
      \code{\link{optimize}}, including the diagonal maximum likelihood
      estimator.}
    \item{\code{etau}}{method-of-moments-like estimator based on
      Kendall's tau for the chosen method.}
  }
}
\author{Marius Hofert}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
  Estimators for Archimedean copulas in high dimensions: A comparison,
  to be submitted.
}
\seealso{
  The more sophisticated estimators \code{\link{emle}} (Maximum Likelihood) and
  \code{\link{emde}} (Minimum Distance). \code{\link{enacopula}}
  (wrapper for different estimators).
}
\examples{
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

system.time(theta.hat.beta <- ebeta(U, cop=cop))
theta.hat.beta$root

system.time(theta.hat.dmle <- edmle(U, cop=cop))
theta.hat.dmle$minimum

system.time(theta.hat.etau <- etau(U, cop=cop, method="tau.mean"))
theta.hat.etau

system.time(theta.hat.etau. <- etau(U, cop=cop, method="theta.mean"))
theta.hat.etau.
}
\keyword{models}
