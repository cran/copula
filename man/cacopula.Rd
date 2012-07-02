\name{cacopula}
\alias{cacopula}
\title{Conditional copula function}
\description{
  Compute the conditional copula function
  \eqn{C(u_d\,|\,u_1,\dots,u_{d-1})}{C(u[d]|u[1],...,u[d-1])}
  of \eqn{u_d}{u[d]} given \eqn{u_1,\dots,u_{d-1}}{u[1],...,u[d-1]}.
}
\usage{
cacopula(u, cop, n.MC=0, log=FALSE)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix; the conditioning is done on the values
    in the first \eqn{d-1} columns.}
  \item{cop}{\code{"\linkS4class{outer_nacopula}"} with specified
    parameters (only Archimedean copulas are currently provided).}
  \item{n.MC}{Monte Carlo sample size.}
  \item{log}{if TRUE the logarithm of the conditional copula function is
    returned.}
}
\value{
  \code{\link{numeric}} vector of length \eqn{n} containing the
  conditional copula function of \eqn{u_d}{u[d]} given
  \eqn{u_1,\dots,u_{d-1}}{u[1],...,u[d-1]}.
}
\author{Marius Hofert}
\note{
  For some (but not all) families, this function also makes sense on the
  boundaries (if the corresponding limits can be computed).
}
\seealso{
  \code{\link{acopula-families}}.
}
\examples{
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
d <- 2
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 1000
U <- rnacopula(n, cop)

U. <- cbind(U[,1], cacopula(U, cop=cop)) # should be ~ U[0,1]^2
plot(U.[,1],U.[,2])
}
\keyword{distribution}
