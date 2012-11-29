\name{K}
\alias{K}
\alias{pK}
\alias{qK}
\alias{dK}
\alias{rK}
\title{Kendall Distribution Function for Archimedean Copulas}
\description{
  The distribution function of the \emph{Kendall distribution} of an Archimedean
  copula is defined as
  \deqn{K(u) = P(C(U_1,U_2,\dots,U_d) \le u),}{K(u) = P(C(U[1],U[2],
    \dots,U[d]) <= u),}
  where \eqn{u \in [0,1]}{u in [0,1]}, and the \eqn{d}-dimensional
  \eqn{(U_1,U_2,\dots,U_d)}{(U[1],U[2],\dots,U[d])} is distributed according
  to the copula \eqn{C}. Note that the random variable
  \eqn{C(U_1,U_2,\dots,U_d)}{C(U[1],U[2],...,U[d])} is known as
  \emph{probability integral transform}. Its distribution function
  \eqn{K} is equal to the identity if \eqn{d = 1}, but is non-trivial for
  \eqn{d \ge 2}{d >= 2}.

  \code{pK}, \code{qK}, \code{dK}, and \code{rK} provide the
  distribution function, quantile function, density, and random number
  generator, respectively, for the Kendall distribution of an
  Archimedean copula.
}
\usage{
% K(u, cop, d, n.MC=0, log=FALSE) # is deprecated
pK(u, cop, d, n.MC=0, log=FALSE)
qK(u, cop, d, n.MC=0,
   method=c("default", "simple", "sort", "discrete", "monoH.FC"),
   u.grid, \dots)
dK(u, cop, d, n.MC=0, log=FALSE)
rK(n, cop, d)
}
\arguments{
  \item{u}{evaluation point(s) (have to be in \eqn{[0,1]}).}
  \item{cop}{\code{\linkS4class{acopula}} with specified parameter, or
    (currently for \code{rK} only) a \code{\linkS4class{outer_nacopula}}.}
  \item{d}{dimension (not used when \code{cop} is an
    \code{\linkS4class{outer_nacopula}}).}
  \item{n.MC}{\code{\link{integer}}, if positive, a Monte Carlo approach
    is applied with sample size equal to \code{n.MC} to evaluate the
    generator derivatives involved; otherwise
    (\code{n.MC=0}) the exact formula is used based on the generator
    derivatives as found by Hofert et al. (2011b).}
  \item{log}{logical indicating if the \code{\link{log}} of \eqn{K}
    should be returned.}
  \item{method}{string for the method to compute the quantile function
    of \eqn{K}.  Currently, one of
    \describe{
      \item{"default"}{The default method. Currently chooses
	\code{method="monoH.FC"} with \code{u.grid=0:128/128}.}
      \item{"simple"}{Straightforward root finding based on
	\code{\link{uniroot}}.}
      \item{"sort"}{Root finding based on \code{\link{uniroot}} but
	first sorting \code{u}.}
      \item{"discrete"}{First, \eqn{K} is evaluated at the given grid
	points \code{u.grid} (which should contain 0 and 1). Based on these
	probabilities, quantiles are computed with \code{\link{findInterval}}.}
      \item{"monoH.FC"}{First, \eqn{K} is evaluated at the given grid points
	\code{u.grid}.  A monotone spline is then used to approximate \eqn{K}.
	Based on this approximation, quantiles are computed with
	\code{\link{uniroot}}.}
    }
  }
  \item{u.grid}{(for \code{method="discrete"}:) the grid on which
  \eqn{K} is evaluated, a \code{\link{numeric}} vector.}
  \item{\dots}{additional arguments passed to \code{\link{uniroot}} (for
    \code{method="default"}, \code{method="simple"}, \code{method="sort"}, and
    \code{method="monoH.FC"}) or \code{\link{findInterval}} (for
    \code{method="discrete"}).}
  \item{n}{sample size for \code{rK}.}
}
\details{
  For a completely monotone Archimedean generator \eqn{\psi}{psi},
  \deqn{K(u)=\sum_{k=0}^{d-1}
    \frac{\psi^{(k)}(\psi^{-1}(u))}{k!} (-\psi^{-1}(u))^k,\ u\in[0,1];}{%
    K(u)=sum(k=0,...,d-1) psi^{(k)}(psi^{-1}(u))/k! (-psi^{-1}(u))^k, u
    in [0,1];}
  see Barbe et al. (1996).  The corresponding density is
  \deqn{\frac{(-1)^d\psi^{(d)}(\psi^{-1}(u))}{(d-1)!}
    (-(\psi^{-1})'(u))(\psi^{-1}(u))^{d-1}}{%
    (-1)^d psi^{(d)}(psi^{-1}(u))(-(psi^{-1})'(u))(psi^{-1}(u))^{d-1}/(d-1)!}
}
\value{The distribution function, quantile function, density, and random
  numbers of a Kendall distribution.}
\author{Marius Hofert}
\references{
  Barbe, P., Genest, C., Ghoudi, K., and R\enc{é}{e}millard, B. (1996),
  On Kendall's Process,
  \emph{Journal of Multivariate Analysis} \bold{58}, 197--229.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{gnacopula}}, \code{\link{htrafo}} or \code{\link{emde}}
  (where \code{\link{K}} is used).
}
\examples{
tau <- 0.5
(theta <- copGumbel@iTau(tau)) # 2
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

## compute Kendall distribution function
u <- seq(0,1, length.out = 255)
Ku    <- pK(u, cop=cop@copula, d=d) # exact
Ku.MC <- pK(u, cop=cop@copula, d=d, n.MC=1000) # via Monte Carlo

## build sample from K
set.seed(1)
n <- 200
W <- rK(n, cop)

## plot empirical distribution function based on W
## and the corresponding theoretical Kendall distribution function
## (exact and via Monte Carlo)
plot(ecdf(W), col="blue", xlim=c(0,1), verticals=TRUE,
     main = expression("Empirical"~ F[n]( C(U) ) ~
                       "and its Kendall distribution"~ K(u)),
     do.points=FALSE, asp=1)
abline(0,1, lty=2); abline(h=0:1, v=0:1, lty=3, col="gray")
lines(u, Ku.MC, col="red")
lines(u, Ku, col="black")
legend(.2,.6, expression(F[n],K(u), K[MC](u)),
       col=c("blue","red","black"), lty=1, bty="n",
       xjust = 0.25, yjust = 0.5)

\donttest{% needs a couple of valuable CRAN seconds :-)
### testing qK
uexpr <- quote( 0:63/63 );  u <- eval(uexpr)
d <- 10
cop <- onacopulaL("Gumbel", list(theta = 2, 1:d))

system.time(Ku1 <- qK(u, cop=cop@copula, d=d, method="simple"))
system.time(Ku2 <- qK(u, cop=cop@copula, d=d, method="sort"))
system.time(Ku3 <- qK(u, cop=cop@copula, d=d, method="discrete",
                      u.grid=0:1e4/1e4))
system.time(Ku4 <- qK(u, cop=cop@copula, d=d, method="monoH.FC",
                      u.grid=0:5e2/5e2))

cols <- adjustcolor(c("gray80", "light blue", "royal blue", "purple3"), 0.6)
matplot(u, cbind(Ku1,Ku2,Ku3,Ku4), type="l", lwd=2*4:1, lty = 1:4, col= cols,
        xlab=substitute(u == U, list(U=uexpr)), ylab=expression({K^{-1}}(u)))
legend("topleft", col=cols, lwd=2*4:1, lty=1:4, bty="n", inset=.03,
       legend= paste0("method= ",
             sQuote(c("simple", "sort", "discrete", "monoH.FC"))))
}
}
\keyword{distribution}
