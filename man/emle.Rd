\name{emle}
\title{Maximum Likelihood Estimators for (Nested) Archimedean Copulas}
\alias{emle}
\alias{.emle}
\description{
  Compute (simulated) maximum likelihood estimators for (nested)
  Archimedean copulas.
}
\usage{
emle(u, cop, n.MC=0, optimizer="optimize", method,
     interval=initOpt(cop@copula@name),
     start=list(theta=initOpt(cop@copula@name, interval=FALSE, u=u)),
     \dots)
.emle(u, cop, n.MC=0,
      interval=initOpt(cop@copula@name), \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula, with \eqn{n} the sample size
    and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated
    (currently only non-nested, that is,
    \ifelse{latex}{Archi-medean}{Archimedean} copulas are admitted).}
  \item{n.MC}{\code{\link{integer}}, if positive, \emph{simulated} maximum
    likelihood estimation (SMLE) is used with sample size equal to
    \code{n.MC}; otherwise (\code{n.MC=0}), MLE.  In SMLE, the \eqn{d}th
    generator derivative and thus the copula density is evaluated via
    (Monte Carlo) simulation, whereas MLE uses the explicit formulas for
    the generator derivatives; see the details below.
  }
  \item{optimizer}{a string or \code{NULL}, indicating the optimizer to
    be used, where \code{NULL} means to use \code{\link{optim}} via the
    standard \R function \code{\link[stats4]{mle}()} from package \pkg{stats4},
    whereas the default, \code{"optimize"} uses \code{\link{optimize}} via
    the \R function \code{\link[bbmle]{mle2}()} from package \pkg{bbmle}.}
  \item{method}{only when \code{optimizer} is \code{NULL} or
    \code{"optim"}, the method to be used for \code{\link{optim}}.}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in
    Hofert et al. (2012).}
  \item{start}{\code{\link{list}} of initial values, passed through.}
  \item{\dots}{additional parameters passed to \code{\link{optimize}}.}
}
\details{
   Exact formulas for the generator derivatives were derived in Hofert
   et al. (2012).  Based on these formulas one can compute the
   (log-)densities of the Archimedean copulas.  Note that for some
   densities, the formulas are numerically highly non-trivial to compute
   and considerable efforts were put in to make the computations
   numerically feasible even in large dimensions (see the source code of
   the Gumbel copula, for example).  Both MLE and SMLE showed good
   performance in the simulation study conducted by Hofert et
   al. (2013) including the challenging 100-dimensional case.
   Alternative estimators (see also \code{\link{enacopula}}) often used
   because of their numerical feasibility, might break down in much
   smaller dimensions.

   Note: SMLE for Clayton currently faces serious numerical issues and is
   due to further research.  This is only interesting from a theoretical point
   of view, since the exact derivatives are known and numerically non-critical
   to evaluate.
}
\value{
  \describe{
    \item{emle}{
      an \R object of class \code{"\link[bbmle:mle2-class]{mle2}"} (and
      thus useful for obtaining confidence intervals) with the
      (simulated) maximum likelihood estimator.}
    \item{.emle}{\code{\link{list}} as returned by
      \code{\link{optimize}()} including the maximum likelihood
      estimator (does not confidence intervals but is typically faster).}
  }
}
%\author{Martin Maechler and Marius Hofert}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2013).
  Archimedean Copulas in High Dimensions: Estimators and Numerical
  Challenges Motivated by Financial Applications.
  \emph{Journal de la Soci\enc{é}{e}t\enc{é}{e} Fran\enc{ç}{c}aise de
    Statistique} \bold{154}(1), 25--63.
}
\seealso{
  \code{\link[bbmle]{mle2}} from package \pkg{bbmle} and
  \code{\link[stats4]{mle}} from \pkg{stats4} on which \code{mle2} is
  modeled. \code{\link{enacopula}} (wrapper for different estimators).
  \code{\link{demo}(opC-demo)} and \code{\link{demo}(GIG-demo)} for
  examples of two-parameter families.
}
\examples{
tau <- 0.25
(theta <- copGumbel@iTau(tau)) # 4/3
d <-  20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n,cop)

## Estimation
system.time(efm <- emle(U, cop))
summary(efm) # using bblme's 'mle2' method

## Profile likelihood plot [using S4 methods from bbmle/stats4] :
pfm <- profile(efm)
ci  <- confint(pfm, level=0.95)
ci
stopifnot(ci[1] <= theta, theta <= ci[2])
plot(pfm)               # |z| against theta, |z| = sqrt(deviance)
plot(pfm, absVal=FALSE, #  z  against theta
     show.points=TRUE) # showing how it's interpolated
## and show the true theta:
abline(v=theta, col="lightgray", lwd=2, lty=2)
axis(1, pos = 0, at=theta, label=quote(theta[0]))

## Plot of the log-likelihood, MLE  and  conf.int.:
logL <- function(x) -efm@minuslogl(x)
       # == -sum(copGumbel@dacopula(U, theta=x, log=TRUE))
logL. <- Vectorize(logL)
I <- c(cop@copula@iTau(0.1), cop@copula@iTau(0.4))
curve(logL., from=I[1], to=I[2], xlab=quote(theta),
      ylab="log-likelihood",
      main="log-likelihood for Gumbel")
abline(v = c(theta, efm@coef), col="magenta", lwd=2, lty=2)
axis(1, at=c(theta, efm@coef), padj = c(-0.5, -0.8), hadj = -0.2,
     col.axis="magenta", label= expression(theta[0], hat(theta)[n]))
abline(v=ci, col="gray30", lwd=2, lty=3)
text(ci[2], extendrange(par("usr")[3:4], f= -.04)[1],
     "95\% conf. int.", col="gray30", adj = -0.1)
}
\keyword{models}
