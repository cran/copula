\name{estim.misc}
\title{Various Estimators for (Nested) Archimedean Copulas}
\alias{ebeta}
\alias{edmle}
\alias{etau}
\description{
  Various Estimators for (Nested) Archimedean Copulas, namely,
  \describe{
    \item{ebeta}{Method-of-moments-like estimator based on (a
      multivariate version of) Blomqvist'sbeta.}
    \item{edmle}{Maximum likelihood estimator based on the diagonal of a
      (nested) Archimedean copula.}
    \item{etau}{Method-of-moments-like estimators based on (bivariate)
      Kendall's tau.}
  }
}
\usage{
ebeta(u, cop, interval = initOpt(cop@copula@name), \dots)
edmle(u, cop, interval = initOpt(cop@copula@name), warn=TRUE, \dots)
 etau(u, cop, method = c("tau.mean", "theta.mean"), warn=TRUE, \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the
    sample size and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated
    (currently only Archimedean copulas are provided).}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in Hofert et
    al. (2013).}
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
	\eqn{[0, 1/3]} and in general if at least one of the computed
	pairwise sample	versions of Kendall's tau is negative.}
    }
  }
  \item{\dots}{additional arguments passed to
    \code{\link{corKendall}} (for \code{etau}, but see \sQuote{Details}),
    to \code{\link{optimize}} (for \code{edmle}), or
    to \code{\link{safeUroot}} (for \code{ebeta}).}
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

  For \code{etau}, \code{\link{corKendall}(u, \dots)} is used and if there
  are no \code{\link{NA}}s in \code{u}, by default (if no additional
  arguments are provided), \code{corKendall()} calls the \eqn{O(n log(n))}
  fast \code{\link[pcaPP]{cor.fk}()} from package \CRANpkg{pcaPP}
  instead of the \eqn{O(n^2)} \code{\link{cor}(*, method="kendall")}.
  Conversely, when \code{u} has \code{NA}s, by default,
  \code{\link{corKendall}(u, \dots)} will use
  \code{\link{cor}(u, method="kendall", use = "pairwise")} such that
  \code{etau(u, *)} will work.
  \cr
  Furthermore, \code{method="tau.mean"} means that the average
  of sample versions of Kendall's tau are computed first and then the
  parameter is determined such that the population version of Kendall's
  tau matches this average (if possible); the \code{method="theta.mean"}
  stands for first computing all pairwise Kendall's tau estimators and
  then returning the mean of these estimators.

  For more details, see Hofert et al. (2013).

  Note that these estimators should be used with care; see the
  performance results in Hofert et al. (2013).  In particular,
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
%\author{Marius Hofert and Martin Maechler}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2013).
  Archimedean Copulas in High Dimensions: Estimators and Numerical
  Challenges Motivated by Financial Applications.
  \emph{Journal de la Soci\enc{é}{e}t\enc{é}{e} Fran\enc{ç}{c}aise de
    Statistique}
  \bold{154}(1), 25--63.
}
\seealso{
  \code{\link{corKendall}()}.

  The more sophisticated estimators \code{\link{emle}} (Maximum Likelihood) and
  \code{\link{emde}} (Minimum Distance). \code{\link{enacopula}}
  (wrapper for different estimators).
}
\examples{
tau <- 0.25
(theta <- copGumbel@iTau(tau)) # 4/3 = 1.333..
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

## etau()  in the case of missing values (NA's)
## ------                 ---------------------
##' @title add Missing Values completely at random
##' @param x  matrix or vector
##' @param prob desired probability ("fraction") of missing values (\code{\link{NA}}s).
##' @return x[] with some (100*prob percent) entries replaced by \code{\link{NA}}s.
addNAs <- function(x, prob) {
    np <- length(x)
    x[sample.int(np, prob*np)] <- NA
    x
}

## UM[] := U[] with 5\% missing values
set.seed(7); UM <- addNAs(U, p = 0.05)
mean(is.na(UM)) # 0.05
## This error if x has NA's was happening for  etau(UM, cop=cop)
## before copula version 0.999-17 (June 2017) :
try(eM <- etau(UM, cop=cop, use = "everything"))
        #  --> Error ... NA/NaN/Inf in foreign function call
## The new default:
eM0 <- etau(UM, cop=cop, use = "pairwise")
eM1 <- etau(UM, cop=cop, use = "complete")
##  use = "complete" is really equivalent to dropping all obs. with with missing values:
stopifnot(all.equal(eM1, etau(na.omit(UM), cop=cop), tol = 1e-15))
## but  use = "pairwise" ---> cor(*, use = "pairwise") is much better:
rbind(etau.U = theta.hat.etau, etau.UM.pairwise = eM0, etau.UM.complete = eM1)
}
\keyword{models}
