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
%\author{Marius Hofert and Martin Maechler}
\references{
  Kemp, A. W. (1981),
  Efficient Generation of Logarithmically Distributed Pseudo-Random Variables,
  \emph{Journal of the Royal Statistical Society: Series C (Applied
    Statistics)} \bold{30}, 3, 249--253.
}
\note{
  In the \pkg{copula} package, the Log(p) distribution is needed only
  for generating Frank copula observations, namely in
  \code{\link{copFrank}@V0()}, where \eqn{p = 1 - exp(-\theta)}, i.e.,
  \code{p = -expm1(-theta)} and \code{Ip = exp(-theta)}.

  For large \eqn{\theta} it would be desirable to pass \code{-theta}
  to \code{rlog()} instead of \code{p}.  This has not yet been implemented.
}
\examples{
## Sample n random variates from a Log(p) distribution and plot a
## "histogram"
n <- 1000
p <- .5
X <- rlog(n, p)
table(X) ## distribution on the integers {1, 2, ..}
## ==> The following plot is more reasonable than a  hist(X, prob=TRUE) :
plot(table(X)/n, type="h", lwd=10, col="gray70")

## case closer to numerical boundary:
lV <- log10(V <- rlog(10000, Ip = 2e-9))# Ip = exp(-theta) <==> theta ~= 20
hV <- hist(lV, plot=FALSE)
dV <- density(lV)
## Plot density and histogram on log scale with nice axis labeling & ticks:
plot(dV, xaxt="n", ylim = c(0, max(hV$density, dV$y)),
     main = "Density of [log-transformed] Log(p),  p=0.999999..")
abline(h=0, lty=3); rug(lV); lines(hV, freq=FALSE, col = "light blue"); lines(dV)
rx <- range(pretty(par("usr")[1:2]))
sx <- outer(1:9, 10^(max(0,rx[1]):rx[2]))
axis(1, at=log10(sx),     labels= FALSE, tcl = -0.3)
axis(1, at=log10(sx[1,]), labels= formatC(sx[1,]), tcl = -0.75)
}
\keyword{distribution}
