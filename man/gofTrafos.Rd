\name{gtrafo}
\alias{gtrafo}
\alias{rtrafo}
\alias{htrafo}
\title{GOF Testing Transformations for Archimedean Copulas}
\description{
  Compute the following goodness-of-fit (GOF) testing transformations,
  \describe{
    \item{Rosenblatt's transformation \code{rtrafo()}}{is also of
      importance outside of GOF computations: The Rosenblatt
      transformation is used for computing conditional copulas and
      for sampling purposes. Currently, \code{rtrafo()} is applicable
      to elliptical and Archimedean copulas.}

    \item{\code{htrafo()}}{the transformation described in Hering and Hofert
      (2014), for Archimedean copulas.}
  }
}
\usage{
rtrafo(u, cop, j.ind = NULL, n.MC=0, inverse=FALSE, log=FALSE)
htrafo(u, cop, include.K=TRUE, n.MC=0, inverse=FALSE,
       method=eval(formals(qK)$method), u.grid, ...)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix with values in \eqn{[0,1]}.  If
    \code{inverse=FALSE} (the default), \code{u} contains
    (pseudo-/copula-)observations from the copula \code{cop} based on
    which the transformation is carried out; consider applying the
    function \code{\link{pobs}()} first in order to obtain \code{u}.  If
    \code{inverse=TRUE}, \code{u} contains \eqn{U[0,1]^d} distributed
    values which are transformed to copula-based (\code{cop}) ones.}
  \item{cop}{a \code{"\linkS4class{Copula}"} with specified parameters
    based on which the transformation is computed.  For \code{htrafo()},
    currently only Archimedean copulas are supported (specified as
    \code{"\linkS4class{outer_nacopula}"}), whereas for \code{rtrafo()},
    \emph{hierarchical} Archimedean and elliptical copulas (see
    \code{\link{ellipCopula}}) are allowed.}
  \item{j.ind}{\code{\link{NULL}} (in which case the Rosenblatt
    transformation is computed (all components)) or an
    \code{\link{integer}} between 2 and \eqn{d}
    indicating the conditional distribution which is to be computed.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{K}} (for \code{htrafo})
    or for approximating the derivatives involved in the Rosenblatt
    transform for Archimedean copulas (for \code{rtrafo}).}
  \item{inverse}{logical indicating whether the inverse of the
    transformation is returned.}
  \item{log}{logical specifying if the logarithm of the transformation,
    i.e., conditional distributions is desired.}
  \item{include.K}{logical indicating whether the last component, involving the
    Kendall distribution function \code{\link{K}}, is used in \code{htrafo}.}
  \item{method}{method to compute \code{qK}.}
  \item{u.grid}{argument of \code{qK} (for \code{method="discrete"}).}
  \item{\dots}{additional arguments passed to \code{qK}.}
}
\details{
  \describe{
    \item{\code{rtrafo}}{Given a \eqn{d}-dimensional random vector
      \eqn{\bm{U}}{U} following an Archimedean copula \eqn{C} with
      generator \eqn{\psi}, the conditional copula of \eqn{U_j=u_j} given
      \eqn{U_1=u_1,\dots,U_{j-1}=u_{j-1}}{U_1=u_1,..., U_{j-1}=u_{j-1}} is
      given by
      \deqn{C(u_j\,|\,u_1,\dots,u_{j-1})=
	\frac{\psi^{(j-1)} \Bigl(\sum_{k=1}^j    \psi^{(-1)}(u_k)\Bigr)}{
	      \psi^{(j-1)} \Bigl(\sum_{k=1}^{j-1}\psi^{(-1)}(u_k)\Bigr)}.}{%
	C(u_j | u_1,...,u_{j-1}) = (psi^{(j-1)}(sum(k=1.. j)  psi^{(-1)}(u_k))) /
                                   (psi^{(j-1)}(sum(k=1..j-1) psi^{(-1)}(u_k))).}
      This formula is either evaluated with the exact derivatives or, if
      \code{n.MC} is positive, via Monte Carlo; see \code{\link{absdPsiMC}}.

      Rosenblatt (1952) showed that
      \eqn{\bm{U}^\prime\sim\mathrm{U}[0,1]^m}{U' ~ U[0,1]^m}, where
      \eqn{U_1^\prime=U_1}{U'_1 = U_1},
      \eqn{U_2^\prime=C(U_2\,|\,U_1)}{U'_2 = C(U_2 | U_1)}, ..., and
      \eqn{U_m^\prime=C(U_m\,|\,U_1,\dots,U_{m-1})}{U'_m = C(U_m | U_1,..., U_{m-1})}.

      \code{rtrafo} applies this transformation row-wise to \code{u}
      (with default \eqn{m=d}) and thus returns an \eqn{n\times m}{n x
	m}-matrix.

      The inverse transformation (\code{inverse=TRUE}) applied to
      \eqn{U[0,1]^d} data is known as \dQuote{conditional distribution
	method} for sampling.

      Note that for the Clayton, the Gauss and the t copula, both the
      conditional copulas and their inverses are known explicitly and
      \code{rtrafo()} utilizes these explicit forms.
    }

    \item{\code{htrafo}}{Given a \eqn{d}-dimensional random vector
      \eqn{\bm{U}}{U} following an Archimedean copula \eqn{C} with
      generator \eqn{\psi}, Hering and Hofert (2014) showed that
      \eqn{\bm{U}^\prime\sim\mathrm{U}[0,1]^d}{U'~U[0,1]^d}, where
      \deqn{U_{j}^\prime=\left(\frac{\sum_{k=1}^{j}\psi^{-1}(U_{k})}{
	  \sum_{k=1}^{j+1}\psi^{-1}(U_{k})}\right)^{j},\ j\in\{1,\dots,d-1\},\
	U_{d}^\prime=K(C(\bm{U})).}{%
	U'_j = ((psi^{-1}(U_1) + ... + psi^{-1}(U_j)) /
	(psi^{-1}(U_1) + ... + psi^{-1}(U_{j+1})))^j,  j in {1,..., d-1},
	U'_d = K(C(U)).}
      \code{htrafo} applies this transformation row-wise to
      \code{u} and thus returns either an \eqn{n\times d}{n x d}- or an
      \eqn{n\times (d-1)}{n x (d-1)}-matrix, depending on whether the last
      component \eqn{U^\prime_d}{U'_d} which involves the (possibly
      numerically challenging) Kendall distribution function \eqn{K} is used
      (\code{include.K=TRUE}) or not (\code{include.K=FALSE}).}
  }
}
\value{
  \code{htrafo()} returns an
  \eqn{n\times d}{n x d}- or \eqn{n\times (d-1)}{n x (d-1)}-matrix
  (depending on whether \code{include.K} is \code{TRUE} or
  \code{FALSE}) containing the transformed input \code{u}.

  \code{rtrafo()} returns an \eqn{n\times d}{n x d}-matrix containing the
  transformed input \code{u}.
}
\author{Marius Hofert, Martin Maechler.}
\references{
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009).
  Goodness-of-fit tests for copulas: A review and a power study.
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--213.

  Rosenblatt, M. (1952).
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics} \bold{23}, 3, 470--472.

  Hering, C. and Hofert, M. (2014).
  Goodness-of-fit tests for Archimedean copulas in high dimensions.
  \emph{Innovations in Quantitative Risk Management}.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{gofCopula}} where both transformations are applied or
  \code{\link{emde}} where \code{htrafo} is applied.
}
\examples{
tau <- 0.5
(theta <- copGumbel@iTau(tau)) # 2
(copG <- onacopulaL("Gumbel", list(theta, 1:5))) # d = 5

set.seed(1)
n <- 1000
x <- rnacopula(n, copG)
x <- qnorm(x) # x now follows a meta-Gumbel model with N(0,1) marginals
u <- pobs(x) # build pseudo-observations

## graphically check if the data comes from a meta-Gumbel model
## with the transformation of Hering and Hofert (2014):
u.h <- htrafo(u, cop=copG) # transform the data
pairs(u.h, gap=0, cex=0.2) # looks good

## with the transformation of Rosenblatt (1952):
u.r <- rtrafo(u, cop=copG) # transform the data
pairs(u.r, gap=0, cex=0.2) # looks good

## what about a meta-Clayton model?
## the parameter is chosen such that Kendall's tau equals (the same) tau
copC <- onacopulaL("Clayton", list(copClayton@iTau(tau), 1:5))

## plot of the transformed data (Hering and Hofert (2014)) to see the
## deviations from uniformity
u.H <- htrafo(u, cop=copC) # transform the data
pairs(u.H, gap=0, cex=0.2) # clearly visible

## plot of the transformed data (Rosenblatt (1952)) to see the
## deviations from uniformity
u.R <- rtrafo(u, cop=copC) # transform the data
pairs(u.R, gap=0, cex=0.2) # clearly visible

## rtrafo() for elliptical:
fN <- fitCopula(normalCopula(dim=ncol(u)), u) # fit a Gauss copula
pairs(rtrafo(u, cop=fN@copula), gap=0, cex=0.2) # visible but not so clearly
if(copula:::doExtras()) {
  f.t <- fitCopula(tCopula(dim=ncol(u)), u)
  tCop <- f.t@copula
} else {
  tCop <- tCopula(param = 0.685, df = 7, dim=ncol(u))
}
u.Rt <- rtrafo(u, cop=tCop) # transform with a fitted t copula
pairs(u.Rt, gap=0, cex=0.2) # *not* clearly visible
\dontshow{%% consistency of  log=FALSE / TRUE:
stopifnot(all.equal(log(u.R),
		    rtrafo(u, cop=copC, log=TRUE), tolerance=1e-14),
	  all.equal(log(rtrafo(u, cop= fN@copula)),
		    rtrafo(u, cop=fN@copula, log=TRUE), tolerance=1e-14),
	  all.equal(log(u.Rt),
		    rtrafo(u, cop=tCop, log=TRUE), tolerance=1e-14),
	  TRUE)
}% only for testing
}
\keyword{transformation}
\keyword{distribution}
\keyword{multivariate}
