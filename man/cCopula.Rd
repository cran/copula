\name{cCopula}
\alias{cCopula}
\alias{cacopula}
\title{Conditional Copula Function}
\description{
  Compute the conditional distribution aka conditional copula function, 
  \eqn{C(u_d\,|\,u_1,\dots,u_{d-1})}{C(u[d]|u[1],...,u[d-1])}
  of \eqn{u_d}{u[d]} given \eqn{u_1,\dots,u_{d-1}}{u[1],...,u[d-1]}.
}
\details{
  \code{cCopula()} can be seen as a special case of the Rosenblatt 
  transform which by default computes all conditional distributions.
  Conseqently, \code{cCopula()} simply calls \code{\link{rtrafo}(...., j.ind=d)}.
}
\usage{
cCopula(u, copula, log=FALSE, n.MC=0)

## Deprecated, for Archimedean copulas only:
cacopula(u, cop, n.MC=0, log=FALSE)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix; the conditioning is done on the values
    in the first \eqn{d-1} columns.}
  \item{copula}{any \code{"\linkS4class{copula}"} with specified
    parameters; currently only Archimedean and elliptical copulas provided.}
  \item{cop}{\code{"\linkS4class{outer_nacopula}"} with specified
    parameters (only Archimedean copulas are currently provided).}
  \item{log}{if TRUE the logarithm of the conditional copula function is
    returned.}
  \item{n.MC}{Monte Carlo sample size; for Archimedean copulas only, used if positive.}
}
\value{
  \code{\link{numeric}} vector of length \eqn{n} containing the
  conditional copula function of \eqn{u_d}{u[d]} given
  \eqn{u_1,\dots,u_{d-1}}{u[1],...,u[d-1]}.
}
\author{Marius Hofert, Yongsheng Wang and Martin Maechler}
\note{
  For some (but not all) families, this function also makes sense on the
  boundaries (if the corresponding limits can be computed).
}
\seealso{
  \code{\link{rtrafo}}; \code{\link{acopula-families}}.
}
\examples{
tau <- 0.5
(theta <- copGumbel@iTau(tau)) # 2
d <- 2
# two ways to specify a Gumbel copula:
cop  <- onacopulaL("Gumbel", list(theta,1:d))
gCop <- gumbelCopula(theta, dim=d) #

n <- 1000
set.seed(1)
U <- rCopula(n, cop)

U. <- cbind(U[,1], cCopula(U, cop=cop)) # should be ~ U[0,1]^2
plot(U.[,1],U.[,2])
\donttest{
 stopifnot(all.equal(cacopula(U, cop),# with deprecation warning
                      cCopula(U, cop)))
}
## more examples:  -->  ?rtrafo % ./gtrafo.Rd
}
\keyword{distribution}
