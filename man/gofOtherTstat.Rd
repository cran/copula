\name{gofOtherTstat}
\alias{gofBTstat}
\title{Various Goodness-of-fit Test Statistics}
\description{
  \code{gofBTstat()} computes supposedly Beta distributed test
  statistics for checking uniformity of \code{u} on the unit sphere.
}
\usage{
gofBTstat(u)
}
\arguments{
  \item{u}{\eqn{(n,d)}-matrix of values whose rows supposedly follow a
    uniform distribution on the unit sphere in \eqn{\mathbf{R}^d}{IR^d}.}
}
\value{
  An \eqn{(n,d-1)}-\code{\link{matrix}} where the \eqn{(i,k)}th entry is
  \deqn{B_{ik}=\frac{\sum_{j=1}^k u_{ij}^2}{\sum_{j=1}^d u_{ij}^2}.
  }{B[ik] = (u[i,1]^2+..+u[i,k]^2) / (u[i,1]^2+..+u[i,d]^2).}
}
%\author{Marius Hofert.}
\references{
  Li, R.-Z., Fang, K.-T., and Zhu, L.-X. (1997).
  Some Q-Q probability plots to test spherical and elliptical symmetry.
  \emph{Journal of Computational and Graphical Statistics} \bold{6}(4), 435--450.
}
\examples{
## generate data on the unit sphere
n <- 360
d <- 5
set.seed(1)
x <- matrix(rnorm(n*d), ncol=d)
U <- x/sqrt(rowSums(x^2))

## compute the test statistics B_k, k in {1,..,d-1}
Bmat <- gofBTstat(U)

## (graphically) check if Bmat[,k] follows a Beta(k/2, (d-k)/2) distribution
qqp <- function(k, Bmat) {
    d <- ncol(Bmat)+1L
    tit <- substitute(plain("Beta")(.A.,.B.)~~ bold("Q-Q Plot"),
                      list(.A. = k/2, .B. = (d-k)/2))
    qqplot2(Bmat[,k], qF=function(p) qbeta(p, shape1=k/2, shape2=(d-k)/2),
            main.args=list(text=tit, side=3, cex=1.3, line=1.1, xpd=NA))
}
qqp(1, Bmat=Bmat) # k=1
qqp(3, Bmat=Bmat) # k=3
}
\keyword{htest}
\keyword{goodness-of-fit}
\keyword{distribution}
\keyword{multivariate}
