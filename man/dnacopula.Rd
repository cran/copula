\name{dnacopula}
\alias{dnacopula}
\title{Density Evaluation for (Nested) Archimedean Copulas}
\description{
  Evaluates the density of an Archimedean copula.
}
\usage{
dnacopula(x, u, log=FALSE, \dots)
}
\arguments{
  \item{x}{an object of class \code{"\linkS4class{outer_nacopula}"}.}
  \item{u}{argument of the copula x.  Note that u can be a matrix in which case
    the density is computed for each row of the matrix and the vector of
    values is returned.}
  \item{log}{logical indicating if the \code{\link{log}} of the density
    should be returned.}
  \item{\dots}{
    optional arguments passed to the copula's \code{dacopula}
    function (slot), such as \code{n.MC} (non-negative integer) for
    possible Monte Carlo evaluation (see \code{dacopula} in
    \code{\linkS4class{acopula}}).}
}
\details{
  If it exists, the density of an Archimedean copula \eqn{C} with
  generator \eqn{\psi}{psi} at \eqn{\bm{u}\in(0,1)^d}{u in (0,1)} is given by
  \deqn{c(\bm{u})=\psi^{(d)}(\psi^{-1}(u_1)+\dots+\psi^{-1}(u_d))
    \prod_{j=1}^d(\psi^{-1}(u_j))^\prime
    = \frac{\psi^{(d)}(\psi^{-1}(u_1)+\dots+\psi^{-1}(u_d))}{
      \prod_{j=1}^d\psi^\prime(\psi^{-1}(u_j))}.
  }{c(u) = psi^{(d)}(psi^{-1}(u_1)+...+psi^{-1}(u_d)) prod(j=1..d)
    (psi^{-1}(u_j))' = psi^{(d)}(psi^{-1}(u_1)+...+psi^{-1}(u_d)) /
     (psi'(psi^{-1}(u_1))*...*psi'(psi^{-1}(u_d))).}
}
\value{
  A \code{\link{numeric}} vector containing the values of the density of the
  Archimedean copula at \code{u}.
}
\author{Marius Hofert, Martin Maechler}
\references{
  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2013).
  Estimators for Archimedean copulas in high dimensions: A comparison.
  submitted.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  For more details about the derivatives of an Archimedean generator,
  see, for example, \code{absdPsi} in class \code{\linkS4class{acopula}}.
}
\examples{
## Construct a twenty-dimensional Gumbel copula with parameter chosen
## such that Kendall's tau of the bivariate margins is 0.25.
theta <- copJoe@iTau(.25)
C20 <- onacopula("J", C(theta, 1:20))

## Evaluate the copula density at the point u = (0.5,...,0.5)
u <- rep(0.5, 20)
dnacopula(C20, u)

## the same with Monte Carlo based on 10000 simulated "frailties"
dnacopula(C20, u, n.MC = 10000)

## Evaluate the exact log-density at several points
u <- matrix(runif(100), ncol=20)
dnacopula(C20, u, log = TRUE)
}
\keyword{distribution}

