\name{pobs}
\alias{pobs}
\title{Pseudo-observations}
\description{
  Compute the pseudo-observations for the given data matrix.
}
\usage{
pobs(x, na.last = "keep",
     ties.method = c("average", "first", "random", "max", "min"))
}
\arguments{
  \item{x}{\eqn{n\times d}{n x d}-matrix of random variates to be converted to
    pseudo-observations.}
  \item{na.last, ties.method}{are passed to \code{\link{rank}}; see there.}
}
\value{matrix of the same dimensions as \code{x} containing the
  pseudo-observations.
}
\details{
  Given \eqn{n} realizations
  \eqn{\bm{x}_i=(x_{i1},\dots,x_{id})^T}{x_i=(x_{i1},...,x_{id})},
  \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} of a random vector \eqn{\bm{X}}{X},
  the pseudo-observations are defined via \eqn{u_{ij}=r_{ij}/(n+1)} for
  \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} and \eqn{j\in\{1,\dots,d\}}{j in
    {1,...,d}}, where \eqn{r_{ij}} denotes the rank of \eqn{x_{ij}} among all
  \eqn{x_{kj}}, \eqn{k\in\{1,\dots,n\}}{k in {1,...,n}}.  The
  pseudo-observations can thus also be computed by component-wise applying the
  empirical distribution functions to the data and scaling the result by
  \eqn{n/(n+1)}.  This asymptotically negligible scaling factor is used to
  force the variates to fall inside the open unit hypercube, for example, to
  avoid problems with density evaluation at the boundaries.
}
\author{Marius Hofert}
\seealso{
  \code{\link{gnacopula}} or \code{\link{enacopula}} (where this function is
  typically applied).
}
\examples{
## Simple definition of the function:
pobs

## Draw from a multivariate normal distribution
d <- 20
P <- matrix(rep(0.75, d*d), ncol=d)
diag(P) <- rep(1, d)
n <- 1000
x <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = P)

## Compute pseudo-observations (should roughly follow a Gaussian
## copula with correlation matrix P)
u <- pobs(x)
plot(u[,1],u[,2])
}


