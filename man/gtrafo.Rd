\name{gtrafo}
\alias{gtrafo}
\alias{htrafo}
\alias{rtrafo}
\title{Goodness-of-Fit Testing Transformations for (Nested) Archimedean Copulas}
\description{
  \code{htrafo()} the transformation described in Hering and Hofert (2011).

  \code{rtrafo()} the transformation of Rosenblatt (1952).
}
\usage{
htrafo(u, cop, include.K=TRUE, n.MC=0, inverse=FALSE,
       method=formals(qK)$method, u.grid, ...)
rtrafo(u, cop, m=d, n.MC=0)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-/copula-)observations (each
    value in \eqn{[0,1]}) from the copula \code{cop} based on which the
    transformation is carried out. Consider applying the function
    \code{\link{pobs}} first in order to obtain \code{u}.}
  \item{cop}{the \code{"\linkS4class{outer_nacopula}"} with specified
    parameters based on which the transformation is computed (currently only
    Archimedean copulas are provided).}
  \item{include.K}{logical indicating whether the last component, involving the
    Kendall distribution function \code{\link{K}}, is used in \code{htrafo}.}
  \item{m}{integer between 2 and \eqn{d} indicating the order up to which the
    conditional distributions are computed (the largest being the
    \eqn{m}th given all with index smaller than \eqn{m}).}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{K}} (for \code{htrafo})
    or for approximating the derivatives involved (for \code{rtrafo}).}
  \item{inverse}{logical indicating whether the inverse of the
    transformation is returned.}
  \item{method}{method to compute \code{qK}.}
  \item{u.grid}{argument of \code{qK} (for \code{method="discrete"}).}
  \item{\dots}{additional arguments passed to \code{qK}.}
}
\details{
  \describe{
    \item{\code{htrafo}}{Given a \eqn{d}-dimensional random vector \eqn{\bm{U}} following an
      Archimedean copula \eqn{C} with generator \eqn{\psi}, Hering and
      Hofert (2011) showed that \eqn{\bm{U}^\prime\sim\mathrm{U}[0,1]^d}{U'~U[0,1]^d}, where
      \deqn{U_{j}^\prime=\left(\frac{\sum_{k=1}^{j}\psi^{-1}(U_{k})}{
	  \sum_{k=1}^{j+1}\psi^{-1}(U_{k})}\right)^{j},\ j\in\{1,\dots,d-1\},\
	U_{d}^\prime=K(C(\bm{U})).}{%
	U'_j=((psi^{-1}(U_1)+...+psi^{-1}(U_j)) /
	(psi^{-1}(U_1)+...+psi^{-1}(U_{j+1})))^j,  j in {1,...,d-1},
	U'_d=K(C(U)).}
      \code{htrafo} applies this transformation row-wise to
      \code{u} and thus returns either an \eqn{n\times d}{n x d}- or an
      \eqn{n\times (d-1)}{n x (d-1)}-matrix, depending on whether the last
      component \eqn{U^\prime_d}{U'_d} which involves the (possibly
      numerically challenging) Kendall distribution function \eqn{K} is used
      (\code{include.K=TRUE}) or not (\code{include.K=FALSE}).}

    \item{\code{rtrafo}}{Given a \eqn{d}-dimensional random vector
      \eqn{\bm{U}} following an Archimedean copula \eqn{C} with
      generator \eqn{\psi}, the conditional copula of \eqn{U_j=u_j} given
      \eqn{U_1=u_1,\dots,U_{j-1}=u_{j-1}}{U_1=u_1,...,U_{j-1}=u_{j-1}} is
      given by
      \deqn{C(u_j\,|\,u_1,\dots,u_{j-1})=
	\frac{\psi^{(j-1)} \Bigl(\sum_{k=1}^j    \psi^{(-1)}(u_k)\Bigr)}{
	      \psi^{(j-1)} \Bigl(\sum_{k=1}^{j-1}\psi^{(-1)}(u_k)\Bigr)}.}{%
	C(u_j | u_1,...,u_{j-1}) = (psi^{(j-1)}(sum(k=1..j)  psi^{(-1)}(u_k))) /
        (psi^{(j-1)}(sum(k=1..j-1) psi^{(-1)}(u_k))).}
      This formula is either evaluated with the exact derivatives or, if
      \code{n.MC} is positive, via Monte Carlo; see \code{\link{psiDabsMC}}.
      Rosenblatt (1952) showed that
      \eqn{\bm{U}^\prime\sim\mathrm{U}[0,1]^m}{U'~U[0,1]^m}, where
      \eqn{U_1^\prime=U_1}{U'_1=U_1},
      \eqn{U_2^\prime=C(U_2\,|\,U_1)}{U'_2=C(U_2|U_1)}, ..., and
      \eqn{U_m^\prime=C(U_m\,|\,U_1,\dots,U_{m-1})}{U'_m=C(U_m|U_1,...,U_{m-1})}.
      \code{rtrafo} applies this transformation row-wise to
      \code{u} (with default \eqn{m=d}) and thus returns an
      \eqn{n\times m}{n x m}-matrix.
    }
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
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009),
  Goodness-of-fit tests for copulas: A review and a power study
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--213.

  Rosenblatt, M. (1952),
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics} \bold{23}, 3, 470--472.

  Hering, C. and Hofert, M. (2011),
  Goodness-of-fit tests for Archimedean copulas in large dimensions,
  submitted.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b),
  Likelihood inference for Archimedean copulas,
  submitted.
}
\seealso{
  \code{\link{gnacopula}} where both transformations are applied or
  \code{\link{emde}} where \code{htrafo} is applied.
}
\examples{
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
(copG <- onacopulaL("Gumbel", list(theta, 1:5))) # d = 5

set.seed(1)
n <- 1000
x <- rnacopula(n, copG)
x <- qnorm(x) # x now follows a meta-Gumbel model with N(0,1) marginals
u <- pobs(x) # build pseudo-observations

## graphically check if the data comes from a meta-Gumbel model
## with the transformation of Hering and Hofert (2011):
u.h <- htrafo(u, cop=copG) # transform the data
pairs(u.h, cex=0.2) # looks good

## with the transformation of Rosenblatt (1952):
u.r <- rtrafo(u, cop=copG) # transform the data
pairs(u.r, cex=0.2) # looks good

## what about a meta-Clayton model?
## the parameter is chosen such that Kendall's tau equals (the same) tau
copC <- onacopulaL("Clayton", list(copClayton@tauInv(tau), 1:5))

## plot of the transformed data (Hering and Hofert (2011)) to see the
## deviations from uniformity
u.prime <- htrafo(u, cop=copC) # transform the data
pairs(u.prime, cex=0.2) # clearly visible

## plot of the transformed data (Rosenblatt (1952)) to see the
## deviations from uniformity
u.prime. <- rtrafo(u, cop=copC) # transform the data
pairs(u.prime., cex=0.2) # clearly visible
}
\keyword{transformation}
\keyword{distribution}
\keyword{multivariate}
