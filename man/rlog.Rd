\name{rlog}
\alias{rlog}
\alias{rlogR}
\title{Sampling Logarithmic Distributions}
\description{
  Generating random variates from a Log(p) distribution with probability
  mass function
  \deqn{p_k=\frac{p^k}{-\log(1-p)k},\ k\in\mathbf{N},
      }{p_k = p^k/(-log(1-p)k), k in IN,}
  where \eqn{p\in(0,1)}{p in (0,1)}.  The implemented algorithm is the
  one named \dQuote{LK} in Kemp (1981).
}
\usage{
   rlog(n, p, Ip = 1 - p)
}
\arguments{
  \item{n}{sample size, that is, length of the resulting vector of random
    variates.}
  \item{p}{parameter in \eqn{(0,1)}.}
  \item{Ip}{\eqn{= 1 - p}, possibly more accurate, e.g, when
    \eqn{p\approx 1}{p ~= 1}.}
}
\value{
  A vector of positive \code{\link{integer}}s of length \code{n} containing the
  generated random variates.
}
\details{
  For documentation and didactical purposes, \code{rlogR} is a pure-\R
  implementation of \code{rlog}.  However, \code{rlogR} is not as fast as
  \code{rlog} (the latter being implemented in C).
}
\author{Marius Hofert, Martin Maechler}
\references{
  Kemp, A. W. (1981),
  Efficient Generation of Logarithmically Distributed Pseudo-Random Variables,
  \emph{Journal of the Royal Statistical Society: Series C (Applied
    Statistics)} \bold{30}, 3, 249--253.
}
\examples{
## Sample n random variates from a Log(p) distribution and plot a
## histogram
n <- 1000
p <- .5
X <- rlog(n, p)
hist(X, prob = TRUE)
}
\keyword{distribution}
