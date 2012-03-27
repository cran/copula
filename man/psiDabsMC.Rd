\name{psiDabsMC}
\alias{psiDabsMC}
\title{Absolute Value of Generator Derivatives via Monte Carlo}
\description{
  Computes the absolute values of the \eqn{d}th generator derivative
  \eqn{\psi^{(d)}}{psi^{(d)}} via Monte Carlo simulation.
}
\usage{
psiDabsMC(t, family, theta, degree = 1, n.MC,
          method = c("log", "direct", "pois.direct", "pois"),
          log = FALSE)
}
\arguments{
  \item{t}{\code{\link{numeric}} vector of evaluation points.}
  \item{family}{Archimedean family (name or object).}
  \item{theta}{parameter value.}
  \item{degree}{order \eqn{d} of the derivative.}
  \item{n.MC}{Monte Carlo sample size.}
  \item{method}{different methods:
    \describe{
      \item{\code{"log"}:}{evaluates the logarithm of the sum involved
	in the Monte Carlo approximation in a numerically stable way;}
      \item{\code{"direct"}:}{directly evaluates the sum;}
      \item{\code{"pois.direct"}:}{interprets the sum in terms of the
	density of a Poisson distribution and evaluates this density directly;}
      \item{\code{"pois"}:}{as for \code{method="pois"} but evaluates
	the logarithm of the Poisson density in a numerically stable way.}
    }
  }
  \item{log}{if TRUE the logarithm of psiDabs is returned.}
}
\details{
  The absolute value of the \eqn{d}th derivative of the Laplace-Stieltjes transform
  \eqn{\psi=\mathcal{LS}[F]}{psi=LS[F]} can be approximated via
  \deqn{(-1)^d\psi^{(d)}(t)=\int_0^\infty
    x^d\exp(-tx)\,dF(x)\approx\frac{1}{N}\sum_{k=1}^NV_k^d\exp(-V_kt),\ t > 0,}{%
    (-1)^d psi^{(d)}(t) = int_0^Inf x^d exp(-tx) dF(x) ~= (1/N) sum(k=1..N)V_k^d exp(-V_k t), t > 0,}
  where \eqn{V_k\sim F,\ k\in\{1,\dots,N\}}{V_k ~ F, k in {1,...,N}}.
  This approximation is used where \eqn{d=}\code{degree} and
  \eqn{N=}\code{n.MC}.  Note that this is comparably fast even if
  \code{t} contains many evaluation points, since the random variates
  \eqn{V_k\sim F,\ k\in\{1,\dots,N\}}{V_k ~ F, k in {1,...,N}} only have
  to be generated once, not depending on \code{t}.
}
\value{
  \code{\link{numeric}} vector of the same length as \code{t} containing
  the absolute values of the generator derivatives.
}
\author{Marius Hofert}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
  Estimators for Archimedean copulas in high dimensions: A comparison,
  to be submitted.
}
\seealso{
  \code{\link{acopula-families}}.
}
\examples{
t <- c(0:100,Inf)
set.seed(1)
psiDabsMC(t, family="Gumbel", theta=2, degree=10, n.MC=10000, log=TRUE)
## Note: The absolute value of the derivative at 0 should be Inf for
## Gumbel, however, it is always finite for the Monte Carlo approximation
}
\keyword{distribution}
