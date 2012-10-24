\name{pnacopula}
\alias{pnacopula}
\title{Evaluation of (Nested) Archimedean Copulas}
\description{
  \code{pnacopula} evaluates a (nested) Archimedean copula (object of class
  \code{\linkS4class{nacopula}}) at the given vector or matrix \code{u}.
}
\usage{
pnacopula(x, u)
}
\arguments{
  \item{x}{(nested) Archimedean copula of dimension \eqn{d}, that is, an
    object of class \code{\linkS4class{nacopula}}, typically from
    \code{\link{onacopula}(..)}.}
  \item{u}{a \code{\link{numeric}} vector of length \eqn{d} or matrix
    with \eqn{d} columns.}
}
\details{
  The value of an Archimedean copula \eqn{C} with generator \eqn{\psi}{psi} at
  \eqn{u} is given by
  \deqn{C(\bm{u})=\psi(\psi^{-1}(u_1)+\dots+\psi^{-1}(u_d)),\ \bm{u}\in[0,1]^d.
  }{    C(u) = psi(psi^{-1}(u_1)+...+psi^{-1}(u_d)), u in [0,1]^d.}
  The value of a nested Archimedean copula is defined similarly.  Note that a
  d-dimensional copula is called \emph{nested Archimedean} if it is an
  Archimedean copula with arguments possibly replaced by other nested
  Archimedean copulas.
}
\value{
  A \code{\link{numeric}} in \eqn{[0,1]} which is the copula evaluated
  at \code{u}. (Currently not parallelized.)
}
\author{Marius Hofert, Martin Maechler.}
\examples{
## Construct a three-dimensional nested Joe copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copJoe@iTau(.2)
theta1 <- copJoe@iTau(.5)
C3 <- onacopula("J", C(theta0, 1, C(theta1, c(2,3))))

## Evaluate this copula at the vector u
u <- c(.7,.8,.6)
pnacopula(C3, u)

## Evaluate this copula at the matrix v
v <- matrix(runif(300), ncol=3)
pnacopula(C3, v)
}
\keyword{multivariate}
\keyword{distribution}
