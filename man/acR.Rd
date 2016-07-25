\name{acR}
\alias{pacR}
\alias{qacR}
\title{Distribution of the Radial Part of an Archimedean Copula}
\description{
  \code{pacR()} computes the distribution function \eqn{F_R} of the radial
  part of an Archimedean copula, given by
  \deqn{F_R(x)=1-\sum_{k=0}^{d-1}
    \frac{(-x)^k\psi^{(k)}(x)}{k!},\ x\in[0,\infty);}{%
    F_R(x)=1-sum(k=0,...,d-1) (-x)^k psi^{(k)}(x)/k!, u in [0,Inf)}
  The formula (in a slightly more general form) is given by
  McNeil and G. \enc{Nešlehová}{Neslehova} (2009).

  \code{qacR()} computes the quantile function of \eqn{F_R}.
}
\usage{
pacR(x, family, theta, d, lower.tail = TRUE, log.p = FALSE, \dots)
qacR(p, family, theta, d, log.p = FALSE, interval,
     tol = .Machine$double.eps^0.25, maxiter = 1000, \dots)
}
\arguments{
  \item{x}{numeric vector of nonnegative evaluation points for \eqn{F_R}.}
  \item{p}{numeric vector %, in \eqn{[0,1]},
    of evaluation points of the  quantile function.}
  \item{family}{Archimedean family.}
  \item{theta}{parameter \eqn{theta}.}
  \item{d}{dimension \eqn{d}.}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE},
    probabilities are \eqn{P[X <= x]} otherwise, \eqn{P[X > x]}.}
  \item{log.p}{\code{\link{logical}}; if \code{TRUE}, probabilities
    \eqn{p} are given as \eqn{\log p}{log(p)}.}
  \item{interval}{root-search interval.}
  \item{tol}{see \code{\link{uniroot}()}.}
  \item{maxiter}{see \code{\link{uniroot}()}.}
  \item{\dots}{additional arguments passed to the procedure for
    computing derivatives.}
}
\value{The distribution function of the radial part evaluated at
  \code{x}, or its inverse, the quantile at \code{p}.}
%\author{Marius Hofert}
\references{
  McNeil, A. J., G. \enc{Nešlehová}{Neslehova}, J. (2009).
  Multivariate Archimedean copulas, \eqn{d}-monotone functions and
  \eqn{l_1}-norm symmetric distributions. \emph{The Annals of Statistics}
  \bold{37}(5b), 3059--3097.
}
\examples{
## setup
family <- "Gumbel"
tau <- 0.5
m <- 256
dmax <- 20
x <- seq(0, 20, length.out=m)

## compute and plot pacR() for various d's
y <- vapply(1:dmax, function(d)
            pacR(x, family=family, theta=iTau(archmCopula(family), tau), d=d),
            rep(NA_real_, m))
plot(x, y[,1], type="l", ylim=c(0,1),
     xlab = expression(italic(x)),
     ylab = substitute(italic(F[R](x))~~"for d=1:"*dm, list(dm=dmax)))
for(k in 2:dmax) lines(x, y[,k])
}
\keyword{distribution}
