\name{pobs}
\alias{pobs}
\title{Pseudo-Observations}
\description{
  Compute the pseudo-observations for the given data matrix.
}
\usage{
pobs(x, na.last = "keep",
     ties.method = eval(formals(rank)$ties.method), lower.tail = TRUE)
}
\arguments{
  \item{x}{\eqn{n\times d}{n x d}-matrix (or \eqn{d}-vector) of random
    variates to be converted to pseudo-observations.}
  \item{na.last}{string passed to \code{\link{rank}}; see there.}
  \item{ties.method}{string specifying how ranks should be computed if
     there are ties in any of the coordinate samples of \code{x}; passed
     to \code{\link{rank}}.}
  \item{lower.tail}{\code{\link{logical}} which, if \code{FALSE},
    returns the pseudo-observations when applying the
    empirical marginal survival functions.}
}
\value{\code{\link{matrix}} (or \code{\link{vector}}) of the same
  dimensions as \code{x} containing the pseudo-observations.
}
\details{
  Given \eqn{n} realizations
  \eqn{\bm{x}_i=(x_{i1},\dots,x_{id})^T}{x_i=(x_{i1},...,x_{id})},
  \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} of a random vector \eqn{\bm{X}}{X},
  the pseudo-observations are defined via \eqn{u_{ij}=r_{ij}/(n+1)} for
  \eqn{i\in\{1,\dots,n\}}{i in {1,...,n}} and \eqn{j\in\{1,\dots,d\}}{j in
    {1,...,d}}, where \eqn{r_{ij}} denotes the rank of \eqn{x_{ij}} among all
  \eqn{x_{kj}}, \eqn{k\in\{1,\dots,n\}}{k in {1,...,n}}.  When there are
  no ties in any of the coordinate samples of \code{x}, the
  pseudo-observations can thus also be computed by component-wise applying the
  marginal empirical distribution functions to the data and scaling the result by
  \eqn{n/(n+1)}. This asymptotically negligible scaling factor is used to
  force the variates to fall inside the open unit hypercube, for example, to
  avoid problems with density evaluation at the boundaries. Note that
  \code{pobs(, lower.tail=FALSE)} simply returns \code{1-pobs()}.
}
%\author{Marius Hofert}
\examples{
## Simple definition of the function:
pobs
%packageDescription("Matrix")% for debugging checks / Matrix
## Draw from a multivariate normal distribution
d <- 10
set.seed(1)
P <- Matrix::nearPD(matrix(pmin(pmax(runif(d*d), 0.3), 0.99), ncol=d))$mat
diag(P) <- rep(1, d)
n <- 500
x <- MASS::mvrnorm(n, mu = rep(0, d), Sigma = P)

## Compute pseudo-observations (should roughly follow a Gauss
## copula with correlation matrix P)
u <- pobs(x)
plot(u[,5],u[,10], xlab=quote(italic(U)[1]), ylab=quote(italic(U)[2]))
\dontshow{if(requireNamespace("zoo")) {% check
 z <- zoo::zoo(x)
 uz <- pobs(z)
 stopifnot( uz == u, identical(dim(uz), dim(u)), inherits(uz, "zoo"))
}}
## All components: pairwise plot
pairs(u, gap=0, pch=".", labels =
      as.expression( lapply(1:d, function(j) bquote(italic(U[.(j)]))) ))
}


