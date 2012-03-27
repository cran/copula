\name{K}
\alias{K}
\title{Kendall Distribution Function of an Archimedean Copula}
\description{
  Compute the \emph{Kendall distribution function} of an Archimedean copula,
  defined as
  \deqn{K(u) = P(C(U_1,U_2,\dots,U_d) \le u),}{K(u) = P(C(U1,U2,\dots,Ud) <= u),}
  where \eqn{u \in [0,1]}{u in [0,1]}, and the \eqn{d}-dimensional
  \eqn{(U_1,U_2,\dots,U_d)}{(U1,U2,\dots,Ud)} is distributed according
  to the copula \eqn{C}.

  Note that the random variable \eqn{C(U_1,U_2,\dots,U_d)}{C(U_1,U_2,...,U_d)}
 is known as \emph{probability integral transform}. Its distribution function
  \eqn{K} is equal to the identity if \eqn{d = 1}, but is non-trivial for
 \eqn{d \ge 2}{d >= 2}.
}
\usage{
K(u, cop, d, n.MC=0)
}
\arguments{
  \item{u}{evaluation point(s) (have to be in \eqn{[0,1]}).}
  \item{cop}{acopula with specified parameter.}
  \item{d}{dimension.}
  \item{n.MC}{\code{\link{integer}}, if positive, a Monte Carlo approach
    is applied with sample size equal to \code{n.MC}; otherwise
    (\code{n.MC=0}) the exact formula is used based on the generator
    derivatives as found by Hofert et al. (2011b).}
}
\details{
  For a completely monotone Archimedean generator \eqn{\psi}{psi},
  \deqn{K(u)=\sum_{k=0}^{d-1}
    \frac{\psi^{(k)}(\psi^{-1}(u))}{k!} (-\psi^{-1}(u))^k,\ u\in[0,1];}{%
    K(u)=sum(k=0,...,d-1) psi^{(k)}(psi^{-1}(u))/k! (-psi^{-1}(u))^k, u in [0,1];}
  see Barbe et al. (1996).
}
\value{Kendall distribution function at \code{u}}
\author{Marius Hofert}
\references{
  Barbe, P., Genest, C., Ghoudi, K., and R\enc{é}{e}millard, B. (1996),
  On Kendall's Process,
  \emph{Journal of Multivariate Analysis}, \bold{58}, 197--229.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b),
  Likelihood inference for Archimedean copulas,
  submitted.
}
\seealso{
  \code{\link{gnacopula}}, \code{\link{htrafo}} or \code{\link{emde}}
  (where \code{\link{K}} is used).
}
\examples{
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

## compute Kendall distribution function
u <- seq(0,1, length.out = 255)
Ku <- K(u, cop=cop@copula, d=d) # exact
Ku.MC <- K(u, cop=cop@copula, d=d, n.MC=1000) # via Monte Carlo

## build sample from K
set.seed(1)
n <- 200
U <- rnacopula(n, cop)
W <- pnacopula(cop, u=U)

## plot empirical distribution function based on W
## and the corresponding theoretical Kendall distribution function
## (exact and via Monte Carlo)
plot(ecdf(W), col="blue", xlim=c(0,1), verticals=TRUE,
     main = expression("Empirical"~ F[n]( C(U) ) ~
                       "and its Kendall distribution"~K(u)),
     do.points=FALSE, asp=1)
abline(0,1, lty=2); abline(h=0:1, v=0:1, lty=3, col="gray")
lines(u, Ku.MC, col="red")
lines(u, Ku, col="black")
legend(.2,.6, expression(F[n],K(u), K[MC](u)),
       col=c("blue","red","black"), lty=1, bty="n",
       xjust = 0.25, yjust = 0.5)
}
\keyword{distribution}
