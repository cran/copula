\name{gofMMDtest}
\alias{gofMMDtest}
\title{Goodness-of-fit Two-Sample Test}
\description{
  \code{gofMMDtest()} computes the two-sample goodness of fit test
  of Schrab et al. (2024), a multiplier bootstrap based on the squared MMD.
}
\usage{
gofMMDtest(x, y, N = 1000, bandwidth2 = 10^c(-4, -3, -2, -3/2, -5/4, -9/8))
}
\arguments{
  \item{x, y}{\eqn{n\times d}{n x d}-matrices of (typically copula) samples
    to be compare. The two matrices must have an equal number
    of rows \eqn{n} and columns \eqn{d}.}
  \item{N}{number of bootstrap replications.}
  \item{bandwidth2}{squared bandwidths of the underlying Gaussian
    mixture kernel based on the average norm (norm scaled by \eqn{1/d}).}
}
\details{
  The test is a multiplier bootstrap test, based on the squared MMD, where
  the underlying kernel is a Gaussian mixture kernel. The scaling of
  the norm in the kernel is done for comparison of test statistics in
  different dimensions.

  Note that his is currently an experimental function, there are
  no known theoretical guarantees for applying this test to
  pseudo-observations, etc.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list with
  components
  \item{method}{a method string identifying the type of test.}
  \item{statistic}{the realized test statistic.}
  \item{statistics.H0}{the bootstrapped test statistics.}
  \item{p.value}{the bootstrapped p-value.}
  \item{data.name}{a string identifying the data used.}
}
\author{Marius Hofert}
\references{
  Schrab, A., Kim, I. Albert, M., Laurent, B., Guedj, B., Gretton, A. (2023)
  MMD aggregated two-sample test. \emph{Journal of Machine Learning Research} \bold{24}, 1--81.
}
\seealso{
  \code{\link{gofT2stat}()} for another two-sample goodness-of-fit
  test statistic.
}
\examples{
n <- 100
d <- 2
set.seed(271)
x <- matrix(runif(n * d), ncol = d)
y <- matrix(runif(n * d), ncol = d)
res <- gofMMDtest(x, y)
res$p.value
}
\keyword{htest}
\keyword{goodness-of-fit}
\keyword{distribution}
\keyword{multivariate}
