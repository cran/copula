\name{enacopula}
\alias{enacopula}
\title{Estimation Procedures for (Nested) Archimedean Copulas}
\description{
  A set of ten different estimators, currently for one-parameter
  Archimedean copulas, of possibly quite high dimensions.
}
\usage{
enacopula(u, cop,
          method = c("mle", "smle", "dmle",
                     "mde.chisq.CvM", "mde.chisq.KS",
                     "mde.gamma.CvM", "mde.gamma.KS",
                     "tau.tau.mean", "tau.theta.mean", "beta"),
          n.MC = if (method == "smle") 10000 else 0,
          interval = initOpt(cop@copula@name),
          xargs = list(), \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula to be estimated, where \eqn{n} denotes
    the sample size and \eqn{d} the dimension.  Consider applying the
    function \code{\link{pobs}} first in order to obtain \code{u}.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated
    (currently only Archimedean copulas are provided).}
  \item{method}{a \code{\link{character}} string specifying the
    estimation method to be used, which has to be one (or a unique
    abbreviation) of
    \describe{
      \item{\code{"mle"}}{maximum likelihood estimator (MLE) computed via \code{\link{.emle}}.}
      \item{\code{"smle"}}{simulated maximum likelihood estimator (SMLE)
	computed with the function \code{\link{.emle}}, where
	\code{n.MC} gives the Monte Carlo sample size.}
      \item{\code{"dmle"}}{MLE based on the diagonal (DMLE); see \code{\link{edmle}}.}
      \item{\code{"mde.chisq.CvM"}}{minimum distance estimator based
	on the chisq distribution and Cram\enc{é}{e}r-von Mises distance; see \code{\link{emde}}.}
      \item{\code{"mde.chisq.KS"}}{minimum distance estimation based on
	the chisq distribution and Kolmogorov-Smirnov distance; see \code{\link{emde}}.}
      \item{\code{"mde.gamma.CvM"}}{minimum distance estimation based on
	the Erlang distribution and Cram\enc{é}{e}r-von Mises distance; see \code{\link{emde}}.}
      \item{\code{"mde.gamma.KS"}}{minimum distance estimation based on
	the Erlang distribution and Kolmogorov-Smirnov distance; see \code{\link{emde}}.}
      \item{\code{"tau.tau.mean"}}{averaged pairwise Kendall's tau estimator}
      \item{\code{"tau.theta.mean"}}{average of pairwise Kendall's tau estimators}
      \item{\code{"beta"}}{multivariate Blomqvist's beta estimator}
    }
  }
  \item{n.MC}{only for \code{method = "smle"}: \code{\link{integer}},
    sample size for simulated maximum likelihood estimation.}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in
    Hofert et al. (2011a).  Used for all methods except
    \code{"tau.tau.mean"} and \code{"tau.theta.mean"}.}
  \item{xargs}{list of additional arguments for the chosen estimation method.}
  \item{\dots}{additional arguments passed to \code{\link{optimize}}.}
}
\details{
  \code{\link{enacopula}} serves as a wrapper for the different
  implemented estimators and provides a uniform framework to utilize
  them.  For more information, see the single estimators as given in the
  section \sQuote{See Also}.

  Note that Hofert, \enc{Mächler}{Maechler}, and McNeil (2011a) compared these
  estimators. Their findings include a rather poor performance and numerically
  challenging problems of some of these estimators. In particular, the
  estimators obtained by \code{method="mde.gamma.CvM"},
  \code{method="mde.gamma.KS"}, \code{method="tau.theta.mean"}, and
  \code{method="beta"} should be used with care (or not at all). Overall, MLE
  performed best (by far).
}
\value{
  the estimated parameter, \eqn{\hat{\theta}}{hat(theta)}, that is, currently a
  number as only one-parameter Archimedean copulas are considered.
}
\author{Marius Hofert, Martin Maechler.}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a).
  Estimators for Archimedean copulas in high dimensions: A comparison,
  to be submitted.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b).
  Likelihood inference for Archimedean copulas,
  submitted.
}
\seealso{
  \code{\link{emle}} which returns an object of \code{"\linkS4class{mle}"}
  providing useful methods not available for other estimators.
  \code{\link{demo}(opC-demo)} and \code{\link{demo}(GIG-demo)} for
  examples of two-parameter families.
  \code{\link{edmle}} for the diagonal maximum likelihood estimator.
  \code{\link{emde}} for the minimum distance estimators.
  \code{\link{etau}} for the estimators based on Kendall's tau.
  \code{\link{ebeta}} for the estimator based on Blomqvist's beta.
}
\examples{
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

(meths <- eval(formals(enacopula)$method)) # "mle", "smle", ...
fun <- function(meth, u, cop, theta){
	run.time <- system.time(val <- enacopula(u, cop=cop, method=meth))
	list(value=val, error=val-theta, utime.ms=1000*run.time[[1]])
}
(res <- sapply(meths, fun, u=U, cop=cop, theta=theta))
}
\keyword{models}
