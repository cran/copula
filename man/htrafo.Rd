\name{htrafo}
\title{GOF Testing Transformation of Hering and Hofert}
\alias{htrafo}
\description{
  The transformation described in Hering and Hofert
  (2014), for Archimedean copulas.
}
\usage{
htrafo(u, copula, include.K = TRUE, n.MC = 0, inverse = FALSE,
       method = eval(formals(qK)$method), u.grid, ...)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix with values in \eqn{[0,1]}.  If
    \code{inverse=FALSE} (the default), \code{u} contains
    (pseudo-/copula-)observations from the copula \code{copula} based on
    which the transformation is carried out; consider applying the
    function \code{\link{pobs}()} first in order to obtain \code{u}.  If
    \code{inverse=TRUE}, \code{u} contains \eqn{U[0,1]^d} distributed
    values which are transformed to copula-based (\code{copula}) ones.}
  \item{copula}{an Archimedean copula specified as
    \code{"\linkS4class{outer_nacopula}"} or \code{"\linkS4class{archmCopula}"}.}
  \item{include.K}{logical indicating whether the last component, involving the
    Kendall distribution function \code{\link{K}}, is used in \code{htrafo()}.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{K}}.}
  \item{inverse}{logical indicating whether the inverse of the
    transformation is returned.}
  \item{method}{method to compute \code{\link{qK}()}.}
  \item{u.grid}{argument of \code{\link{qK}()} (for \code{method="discrete"}).}
  \item{\dots}{additional arguments passed to \code{\link{qK}()} if
    \code{inverse = TRUE}.}
}
\details{
  Given a \eqn{d}-dimensional random vector
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
  (\code{include.K=TRUE}) or not (\code{include.K=FALSE}).
}
\value{
  \code{htrafo()} returns an
  \eqn{n\times d}{n x d}- or \eqn{n\times (d-1)}{n x (d-1)}-matrix
  (depending on whether \code{include.K} is \code{TRUE} or
  \code{FALSE}) containing the transformed input \code{u}.
}
%\author{Marius Hofert and Martin Maechler}
\references{
  Hering, C. and Hofert, M. (2014).
  Goodness-of-fit tests for Archimedean copulas in high dimensions.
  \emph{Innovations in Quantitative Risk Management}.
}
\examples{
## Sample and build pseudo-observations (what we normally have available)
## of a Clayton copula
tau <- 0.5
theta <- iTau(claytonCopula(), tau = tau)
d <- 5
cc <- claytonCopula(theta, dim = d)
set.seed(271)
n <- 1000
U <- rCopula(n, copula = cc)
X <- qnorm(U) # X now follows a meta-Gumbel model with N(0,1) marginals
U <- pobs(X) # build pseudo-observations

## Graphically check if the data comes from a meta-Clayton model
## with the transformation of Hering and Hofert (2014):
U.H <- htrafo(U, copula = cc) # transform the data
splom2(U.H, cex = 0.2) # looks good

## The same for an 'outer_nacopula' object
cc. <- onacopulaL("Clayton", list(theta, 1:d))
U.H. <- htrafo(U, copula = cc.)
splom2(U.H., cex = 0.2) # looks good

## What about a meta-Gumbel model?
## The parameter is chosen such that Kendall's tau equals (the same) tau
gc <- gumbelCopula(iTau(gumbelCopula(), tau = tau), dim = d)

## Plot of the transformed data (Hering and Hofert (2014)) to see the
## deviations from uniformity
U.H.. <- htrafo(U, copula = gc)
splom2(U.H.., cex = 0.2) # deviations visible
}
\keyword{transformation}
\keyword{distribution}
\keyword{multivariate}
