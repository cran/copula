\name{cCopula}
\alias{cacopula}
\alias{rtrafo}
\alias{cCopula}
\title{Conditional Copula Function and Its Inverse}
\description{
  Compute the conditional distribution aka conditional copula function,
  \eqn{C(u_d\,|\,u_1,\dots, u_{d-1})}{C(u[d] | u[1],..., u[d-1])}
  of \eqn{u_d}{u[d]} given \eqn{u_1,\dots,u_{d-1}}{u[1],..., u[d-1]}.
}
\usage{
cCopula(u, copula, indices = 1:dim(copula), inverse = FALSE, log =
        FALSE, ...)

## Deprecated (use cCopula() instead):
rtrafo(u, copula, indices = 1:dim(copula), inverse = FALSE, log = FALSE)
cacopula(u, cop, n.MC = 0, log = FALSE)
}
\arguments{
  \item{u}{A data \code{\link{matrix}} in \eqn{[0,1]^(n, d)} of
    \eqn{U(0,1)^d} samples if \code{inverse = FALSE}
    and (pseudo-/copula-)observations if \code{inverse = TRUE}.}
  \item{copula, cop}{An object of class \code{"\linkS4class{Copula}"} with specified
    parameters; currently only Archimedean and elliptical copulas are
    provided.}
  \item{indices}{A vector of indices \eqn{j} (in \eqn{\{1,\dots,d\}}{{1,..,d}}
    (\eqn{d =} copula dimension); unique; sorted in increasing order) for which
    \eqn{C_{j|1,\dots,j-1}(u_j\,|\,u_1,\dots,u_{j-1})}{C_{j|1,\dots,j-1}(u[j]|u[1],...,u[j-1])} (or, if
    \code{inverse = TRUE},
    \eqn{C^-_{j|1,\dots,j-1}(u_j\,|\,u_1,\dots,u_{j-1})}{C^-_{j|1,\dots,j-1}(u[j]|u[1],...,u[j-1])})
    is computed.}
  \item{inverse}{A \code{\link{logical}} indicating whether the inverse
    \eqn{C^-_{j|1,\dots,j-1}(u_j\,|\,u_1,\dots,u_{j-1})}{C^-_{j|1,\dots,j-1}(u[j]|u[1],...,u[j-1])}
    is returned.}
  \item{n.MC}{Monte Carlo sample size; for Archimedean copulas only,
    used if positive.}
  \item{log}{A \code{\link{logical}} indicating whether the logarithmic
    values are returned.}
  \item{\dots}{Additional arguments (currently only used if
    \code{inverse = TRUE} in which case they are passed on to the
    underlying \code{\link{uniroot}()}).}
}
\value{
  An \eqn{(n, k)}-\code{\link{matrix}} where \eqn{k} denotes the length
  of \code{indices}. This matrix contains the
  conditional copula function values
  \eqn{C_{j|1,\dots,j-1}(u_j\,|\,u_1,\dots,u_{j-1})}{C_{j|1,\dots,j-1}(u[j]|u[1],...,u[j-1])}
  or, if \code{inverse = TRUE}, their inverses
  \eqn{C^-_{j|1,\dots,j-1}(u_j\,|\,u_1,\dots,u_{j-1})}{C^-_{j|1,\dots,j-1}(u[j]|u[1],...,u[j-1])}
  for all \eqn{j} in \code{indices}.
}
\details{
  By default, \code{cCopula()} computes the Rosenblatt
  transform; see Rosenblatt (1952). The involved high-order derivatives
  for Archimedean copulas were derived in Hofert et al. (2012).

  \emph{Sampling}, that is, random number generation,
  can be achieved by using \code{inverse=TRUE}. In this case,
  the inverse Rosenblatt transformation is used, which, for sampling
  purposes, is also known as \emph{conditional distribution method}.
  Note that, for Archimedean copulas not being Clayton, this can be slow
  as it involves numerical root finding in each (but the first) component.
}
%\author{Marius Hofert, Yongsheng Wang, and Martin Maechler}
\note{
  For some (but not all) families, this function also makes sense on the
  boundaries (if the corresponding limits can be computed).
}
\references{
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009).
  Goodness-of-fit tests for copulas: A review and a power study.
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--213.

  Rosenblatt, M. (1952).
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics} \bold{23}, 3, 470--472.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{htrafo}}; \code{\link{acopula-families}}.
}
\examples{
## 1) Sampling a conditional Clayton copula given u_1

## Define the copula
tau <- 0.5
theta <- iTau(claytonCopula(), tau = tau)
d <- 2
cc <- claytonCopula(theta, dim = d)
n <- 1000
set.seed(271)

## A small u_1
u1 <- 0.05
U <- cCopula(cbind(u1, runif(n)), copula = cc, inverse = TRUE)
plot(U[,2], ylab = expression(U[2]))

## A large u_1
u1 <- 0.95
U <- cCopula(cbind(u1, runif(n)), copula = cc, inverse = TRUE)
plot(U[,2], ylab = expression(U[2]))


## 2) Sample via conditional distribution method and then apply the
##    Rosenblatt transform
##    Note: We choose the numerically more involved (and thus slower)
##          Gumbel case here

## Define the copula
tau <- 0.5
theta <- iTau(gumbelCopula(), tau = tau)
d <- 5
gc <- gumbelCopula(theta, dim = d)
n <- 200
set.seed(271)
U. <- matrix(runif(n*d), ncol = d) # U(0,1)^d

\donttest{
## Transform to Gumbel sample via conditional distribution method
U <- cCopula(U., copula = gc, inverse = TRUE) # slow for ACs except Clayton
splom2(U) # scatter-plot matrix copula sample

## Rosenblatt transform back to U(0,1)^d (as a check)
U. <- cCopula(U, copula = gc)
splom2(U.) # U(0,1)^d again
}% dont

## 3) cCopula() for elliptical copulas

tau <- 0.5
theta <- iTau(claytonCopula(), tau = tau)
d <- 5
cc <- claytonCopula(theta, dim = d)
set.seed(271)
n <- 1000
U <- rCopula(n, copula = cc)
X <- qnorm(U) # X now follows a meta-Clayton model with N(0,1) marginals
U <- pobs(X) # build pseudo-observations
fN <- fitCopula(normalCopula(dim = d), data = U) # fit a Gauss copula
U.RN <- cCopula(U, copula = fN@copula)
splom2(U.RN, cex = 0.2) # visible but not so clearly
f.t <- fitCopula(tCopula(dim = d), U)
tc <- f.t@copula
U.Rt <- cCopula(U, copula = tc) # transform with a fitted t copula
splom2(U.Rt, cex = 0.2) # still visible but not so clear
tol <- 1e-14
stopifnot(
    all.equal(log(U.RN),
              cCopula(U, copula = fN@copula, log = TRUE), tolerance = tol),
    all.equal(log(U.Rt),
              cCopula(U, copula = tc, log = TRUE), tolerance = tol),
    TRUE)
}
\keyword{distribution}
