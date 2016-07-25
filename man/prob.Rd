\name{prob}
\alias{prob}
\alias{prob-methods}
\alias{prob,Copula-method}
\title{Computing Probabilities of Hypercubes}
\description{
  Compute probabilities of a \eqn{d-}dimensional random vector \eqn{U}
  distributed according to a given copula \code{x} to
  fall in a hypercube \eqn{(l,u]}, where \eqn{l} and \eqn{u} denote the
  lower and upper corners of the hypercube, respectively.
}
\usage{
prob(x, l, u)
}
\arguments{
  \item{x}{copula of dimension \eqn{d}, that is, an object
    inheriting from \code{\linkS4class{Copula}}.}
  \item{l, u}{\eqn{d}-dimensional, \code{\link{numeric}}, lower and
    upper hypercube boundaries, respectively, satisfying
    \eqn{0 \le l_i \le u_i \le 1}{0 <= l_i <= u_i <= 1},
    for \eqn{i\in{1,\dots,d}}{i in {1,...,d}}.}
}
\value{A \code{\link{numeric}} in \eqn{[0,1]} which is the probability
  \eqn{P(l_i< U_i \le u_i)}{P(l[i] < U[i] <= u[i])}.
}
%\author{Marius Hofert and Martin Maechler}
\seealso{
  \code{\link{pCopula}(.)}.
}
\examples{
## Construct a three-dimensional nested Joe copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copJoe@iTau(.2)
theta1 <- copJoe@iTau(.5)
C3 <- onacopula("J", C(theta0, 1, C(theta1, c(2,3))))

## Compute the probability of a random vector distributed according to
## this copula to fall inside the cube with lower point l and upper
## point u.
l <- c(.7,.8,.6)
u <- c(1,1,1)
prob(C3, l, u)

## ditto for a bivariate normal copula with rho = 0.8 :
prob(normalCopula(0.8), c(.2,.4), c(.3,.6))
}
\keyword{distribution}
\keyword{methods}

