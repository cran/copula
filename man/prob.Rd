\name{prob}
\alias{prob}
\alias{prob-methods}
\alias{prob,outer_nacopula-method}
\title{Computing Probabilities of "nacopula" Objects to Fall in Hypercubes}
\description{
  Compute probabilities of a \eqn{d-}dimensional random vector \eqn{U}
  distributed according to the nested Archimedean copula \code{x} to
  fall in a hypercube \eqn{(l,u]}, where \eqn{l} and \eqn{u} denote the
  lower and upper corners of the hypercube, respectively.
}
\usage{
prob(x, l, u)
}
\arguments{
  \item{x}{(nested) Archimedean copula of dimension \eqn{d}, that is, an object
    of class \code{\linkS4class{nacopula}}, typically constructed with
    \code{\link{onacopula}(..)}.}
  \item{l, u}{\eqn{d}-dimensional, \code{\link{numeric}}, lower and
    upper hypercube boundaries, respectively, satisfying
    \eqn{0 \le l_i \le u_i \le 1}{0 <= l_i <= u_i <= 1},
    for \eqn{i\in{1,\dots,d}}{i in {1,...,d}}.}
}
\value{A \code{\link{numeric}} in \eqn{[0,1]} which is the probability
  \eqn{P(l_i< U_i \le u_i)}{P(l[i] < U[i] <= u[i])}.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
  \code{\link{onacopula}()} for the construction of nested Archimedean
  copulas, the class \code{\linkS4class{nacopula}}, and the evaluation
  of nested Archimedean copulas via \code{\link{pnacopula}(.)}.
}
\examples{
## Construct a three-dimensional nested Joe copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copJoe@tauInv(.2)
theta1 <- copJoe@tauInv(.5)
C3 <- onacopula("J", C(theta0, 1, C(theta1, c(2,3))))

## Compute the probability of a random vector distributed according to
## this copula to fall inside the cube with lower point l and upper
## point u.
l <- c(.7,.8,.6)
u <- c(1,1,1)
prob(C3, l, u)
}
\keyword{distribution}
\keyword{methods}

