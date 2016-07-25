\name{initOpt}
\alias{initOpt}
\title{Initial Interval or Value for Parameter Estimation of Archimedean
  Copulas}
\description{
  Compute an initial interval or initial value for optimization/estimation
  routines (only a heuristic; if this fails, choose your own interval or value).
}
\usage{
initOpt(family, tau.range=NULL, interval = TRUE, u,
        method = c("tau.Gumbel", "tau.mean"), warn = TRUE, \dots)
}
\arguments{
  \item{family}{Archimedean family to find an initial interval for.}
  \item{tau.range}{numeric vector containing lower and upper admissible
    Kendall's tau, or \code{NULL} which choses family-specific defaults,
    see the function definition.}
  \item{interval}{\code{\link{logical}} indicating whether an initial interval
    (the default) or an initial value should be returned.}
  \item{u}{matrix of realizations following the copula family specified
    by \code{family}.  Note that \code{u} can be omitted if
    \code{interval=TRUE}.}
  \item{method}{a \code{\link{character}} string specifying the method to be
    used to compute an estimate of Kendall's tau.  This has to be one (or a
    unique abbreviation) of
    \describe{
      \item{\code{"tau.Gumbel"}}{an estimator based on the diagonal
	maximum-likelihood estimator for Gumbel is used.}
      \item{\code{"tau.mean"}}{an estimator based on the mean of pairwise sample
	versions of Kendall's tau is applied.}
    }
  }
  \item{warn}{logical indicating if warnings are printed for
    \code{method="tau.Gumbel"} when the diagonal maximum-likelihood
    estimator is smaller than \eqn{1}.}
  \item{\dots}{additional arguments passed to \code{cor()} when
    \code{method="tau.mean"}.  Note that otherwise (no additional arg.),
    the much faster \code{\link[pcaPP]{cor.fk}()} from package
    \CRANpkg{pcaPP} is used.}
}
\details{
  For \code{method="tau.mean"} and \code{interval=FALSE}, the mean of pairwise
  sample versions of Kendall's tau is computed as an estimator of the Kendall's
  tau of the Archimedean copula family provided.  This can be slow, especially
  if the dimension is large.  Method \code{method="tau.Gumbel"} (the default)
  uses the explicit and thus very fast diagonal maximum-likelihood estimator for
  Gumbel's family to find initial values.  Given this estimator
  \eqn{\hat{\theta}^\mathrm{G}}{hat(theta)^G}, the corresponding Kendall's tau
  is \eqn{\tau^\mathrm{G}(\hat{\theta}^\mathrm{G})}{tau^G(hat(theta)^G)} where
  \eqn{\tau^\mathrm{G}(\theta)=(\theta-1)/\theta}{tau^G(theta)=(theta-1)/theta}
  denotes Kendall's tau for Gumbel.  This provides an estimator of Kendall's tau
  which is typically much faster to evaluate than, pairwise Kendall's taus.
  Given the estimated \sQuote{amount of concordance} based on Kendall's tau, one
  can obtain an initial value for the provided family by applying
  \eqn{\tau^{-1}}{tau^(-1)}, that is, the inverse of Kendall's tau of the
  family for which the initial value is to be computed.  Note that if the
  estimated Kendall's tau does not lie in the range of Kendall's tau as provided
  by the bivariate vector \code{tau.range}, the point in \code{tau.range}
  closest to the estimated Kendall's tau is chosen.

  The default (\code{interval=TRUE}) returns a reasonably large initial
  interval; see the default of \code{tau.range} in the definition of
  \code{initOpt} for the chosen values (in terms of Kendall's tau).  These
  default values cover a large range of concordance.  If this interval is
  (still) too small, one can adjust it by providing \code{tau.range}. If it is
  too large, a \sQuote{distance to concordance} can be used to determine
  parameter values such that the corresponding Kendall's taus share a certain
  distance to the initial value.  For more details, see Hofert et al. (2012).
  Finally, let us note that for the case \code{interval=TRUE}, \code{u} is not
  required.
}
\value{initial interval which can be used for optimization (for example, for
  \code{\link{emle}}).
}
%\author{Marius Hofert}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{enacopula}}, \code{\link{emle}}, \code{\link{edmle}},
  \code{\link{emde}}, and \code{\link{ebeta}} (where \code{initOpt} is
  applied to find initial intervals).
}
\examples{
## Definition of the function:
initOpt

## Generate some data:
tau <- 0.25
(theta <- copGumbel@iTau(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

## Initial interval:
initOpt("Gumbel") # contains theta

## Initial values:
initOpt("Gumbel", interval=FALSE, u=U) # 1.3195
initOpt("Gumbel", interval=FALSE, u=U, method="tau.mean") # 1.2844
}
