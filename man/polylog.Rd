\name{polylog}
\title{Polylogarithm \eqn{\mathrm{Li_s(z)}}{Li_s(z)} and Debye Functions}
\alias{polylog}
\concept{Polylogarithm}
\concept{Dilogarithm}
\concept{dilog}
\alias{debye1}
\alias{debye2}
\description{
  Compute the polylogarithm function \eqn{\mathrm{Li_s(z)}}{Li_s(z)},
  initially defined as the power series,
  \deqn{\mathrm{Li}_s(z) = \sum_{k=1}^\infty {z^k \over k^s},}{%
    Li_s(z) = sum(k=1..Inf; z^k / k^s),}
  for \eqn{|z| < 1}, and then more generally (by analytic continuation) as
  \deqn{\mathrm{Li}_1(z) = -\log(1-z),}{Li_1(z) = -log(1-z),}
    and
  \deqn{\mathrm{Li}_{s+1}(z) = \int_0^z \frac{\mathrm{Li}_s(t)}{t}\,dt.}{%
    Li_{s+1}(z) = Int[0..z] (Li_s(t) / t) dt.}

  Currently, mainly the case of negative integer \eqn{s} is well supported,
  as that is used for some of the Archimedean copula densities.

  For \eqn{s = 2}, \eqn{\mathrm{Li_2(z)}}{Li_2(z)} is also called
  \sQuote{dilogarithm} or \dQuote{Spence's function}.  The
  \code{"default"} method uses the \code{\link[gsl]{dilog}} or
  \code{\link[gsl]{complex_dilog}} function from package \CRANpkg{gsl},
  respectively when \eqn{s = 2}.

  Also compute the Debye_n functions, for \eqn{n=1} and \eqn{n=2}, in a
  slightly more general manner than the \pkg{gsl} package functions
  \code{\link[gsl]{debye_1}} and \code{debye_2} (which cannot deal with
  non-finite \code{x}.)
}
\usage{
polylog(z, s,
        method = c("default", "sum", "negI-s-Stirling",
                   "negI-s-Eulerian", "negI-s-asymp-w"),
        logarithm = FALSE, is.log.z = FALSE, is.logmlog = FALSE,
        asymp.w.order = 0, n.sum)

debye1(x)
debye2(x)
}
\arguments{
  \item{z}{numeric or complex vector}
  \item{s}{complex number; current implementation is aimed at
    \eqn{s \in \{0,-1,\dots\}}{s in (0,-1,...)}}
  \item{method}{a string specifying the algorithm to be used.}
  \item{logarithm}{logical specified to return log(Li.(.)) instead of Li.(.)}
  \item{is.log.z}{logical; if TRUE, the specified \code{z} argument is
    really \eqn{w = \log(z)}{w = log(z)};
    that is, we compute \eqn{\mathrm{Li}_s(\exp(w))}{Li_s(exp(w))}, and
    we typically have \eqn{w < 0}, or equivalently, \eqn{z < 1}.}
  \item{is.logmlog}{logical; if TRUE, the specified argument \code{z} is
    \eqn{lw = \log(-w) = \log(-\log(z))}{lw = log(-w) = log(-log(z))}
    (where as above, \eqn{w = \log(z)}{w = log(z)}).}
  \item{asymp.w.order}{currently only default is implemented.}
  \item{n.sum}{for \code{method="sum"} only: the number of terms used.}

  \item{x}{numeric vector, may contain \code{Inf}, \code{\link{NA}}, and
    negative values.}
}
\details{
  Almost entirely taken from
  \url{https://en.wikipedia.org/wiki/Polylogarithm}:

For integer values of the polylogarithm order, the following
explicit expressions are obtained by repeated application of
\eqn{z \frac{\partial}{\partial z}}{z * d/dz} to
\eqn{\mathrm{Li}_1(z)}{Li_1(z)}:
%---
\deqn{
  \mathrm{Li}_{1}(z) = -\log(1-z), \ \
  \mathrm{Li}_{0}(z) = {z \over 1-z}, \ \
  \mathrm{Li}_{-1}(z) = {z \over (1-z)^2}, \ \
  \mathrm{Li}_{-2}(z) = {z \,(1+z) \over (1-z)^3},
}{%
  Li_1(z) = -log(1-z),
  Li_0(z) =  z / (1-z),
  Li_{-1}(z) = z / (1-z)^2,
  Li_{-2}(z) = z (1+z) / (1-z)^3,
}
\eqn{\mathrm{Li}_{-3}(z) = {z \,(1+4z+z^2) \over (1-z)^4}}{%
  Li_{-3}(z) = z (1+4z+z^2) / (1-z)^4}, etc.

Accordingly, the polylogarithm reduces to a ratio of polynomials in
z, and is therefore a rational function of z, for all nonpositive
integer orders.  The general case may be expressed as a finite sum:
%---

\deqn{\mathrm{Li}_{-n}(z) =
  \left(z \,{\partial \over \partial z} \right)^n \frac{z}{1-z}
  = \sum_{k=0}^n k! \,S(n+1,k+1) \left({z \over {1-z}} \right)^{k+1}
  \ \ (n=0,1,2,\ldots),}{%
  Li_{-n}(z) = ( z d/dz )^n  z/(1-z)
             = sum(k=0..n ; k! S(n+1,k+1) (z /(1-z))^(k+1)),   (n=0,1,2,...),}
where \eqn{S(n,k)} are the Stirling numbers of the second kind.

Equivalent formulae applicable to negative integer orders are
(Wood 1992, § 6) ...
\deqn{\mathrm{Li}_{-n}(z) = {1 \over (1-z)^{n+1}} \sum_{k=0}^{n-1}
  			\left\langle {n \atop k} \right\rangle z^{n-k} =
			\frac{z \sum_{k=0}^{n-1} \left\langle {n \atop k} \right
			  \rangle z^k}{(1-z)^{n+1}},
      \qquad (n=1,2,3,\ldots) ~, }{%
  Li_{-n}(z) = 1/((1-z)^(n+1)) sum(k=0..(n-1);  < n \ k >  z^(n-k)) =
  = (z \sum_{k=0}^{n-1} < n \ k >  z^k) / ((1-z)^(n+1)),
  (n=1,2,3,..),}
where \eqn{\left\langle {n \atop k} \right\rangle}{< n \ k >}  are the
Eulerian numbers; see also \code{\link{Eulerian}}.

% All roots of Li_{-n}(z) are distinct and real; they include z = 0.

% Duplication formula:  2^{1-s} Li_s(z^2) = Li_s(z) + Li_s(-z).
}
\value{
  numeric/complex vector as \code{z}, or \code{x}, respectively.
}
%\author{Martin Maechler}
\references{
  Wikipedia (2011) \emph{Polylogarithm},
  \url{https://en.wikipedia.org/wiki/Polylogarithm}.

  Wood, D. C. (June 1992).
  The Computation of Polylogarithms.  Technical Report 15-92.
  Canterbury, UK: University of Kent Computing Laboratory.
  \url{https://www.cs.kent.ac.uk/pubs/1992/110/}.% Retrieved 2020-11-16.

  Apostol, T. M. (2010), \emph{"Polylogarithm"}, in the
  NIST Handbook of Mathematical Functions, \url{https://dlmf.nist.gov/25.12}

  Lewin, L. (1981).
  \emph{Polylogarithms and Associated Functions}.
  New York: North-Holland. ISBN 0-444-00550-1.

  For Debye functions:  Levin (1981) above, and \cr
  Wikipedia (2014) \emph{Debye function},
  \url{https://en.wikipedia.org/wiki/Debye_function}.
}
\seealso{
  The polylogarithm is used in MLE for some Archimedean copulas; see
  \code{\link{emle}};

  The Debye functions are used for \code{\link{tau}} or
  \code{\link{rho}} computations of the Frank copula.
}
\examples{
## The dilogarithm,  polylog(z, s = 2) = Li_2(.) -- mathmatically defined on C \ [1, Inf)
## so x -> 1 is a limit case:
polylog(z = 1, s = 2)
## in the limit, should be equal to
pi^2 / 6

## Default method uses  GSL's dilog():
rLi2 <- curve(polylog(x, 2), -5, 1, n= 1+ 6*64, col=2, lwd=2)
abline(c(0,1), h=0,v=0:1, lty=3, col="gray40")
## "sum" method gives the same for |z| < 1 and large number of terms:
ii <- which(abs(rLi2$x) < 1)
stopifnot(all.equal(rLi2$y[ii],
            polylog(rLi2$x[ii], 2, "sum", n.sum = 1e5),
          tolerance = 1e-15))


z1 <- c(0.95, 0.99, 0.995, 0.999, 0.9999)
L   <- polylog(         z1,  s=-3,method="negI-s-Euler") # close to Inf
LL  <- polylog(     log(z1), s=-3,method="negI-s-Euler",is.log.z=TRUE)
LLL <- polylog(log(-log(z1)),s=-3,method="negI-s-Euler",is.logmlog=TRUE)
all.equal(L, LL)
all.equal(L, LLL)

p.Li <- function(s.set, from = -2.6, to = 1/4, ylim = c(-1, 0.5),
                 colors = c("orange","brown", palette()), n = 201, ...)
{
    s.set <- sort(s.set, decreasing = TRUE)
    s <- s.set[1] # <_ for auto-ylab
    curve(polylog(x, s, method="negI-s-Stirling"), from, to,
          col=colors[1], ylim=ylim, n=n, ...)
    abline(h=0,v=0, col="gray")
    for(is in seq_along(s.set)[-1])
        curve(polylog(x, s=s.set[is], method="negI-s-Stirling"),
              add=TRUE, col = colors[is], n=n)
    s <- rev(s.set)
    legend("bottomright",paste("s =",s), col=colors[2-s], lty=1, bty="n")
}

## yellow is unbearable (on white):
palette(local({p <- palette(); p[p=="yellow"] <- "goldenrod"; p}))

## Wikipedia page plot (+/-):
p.Li(1:-3, ylim= c(-.8, 0.6), colors = c(2:4,6:7))

## and a bit more:
p.Li(1:-5)

## For the range we need it:
ccol <- c(NA,NA, rep(palette(),10))
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(0, 10))
## log-y scale:
p.Li(-1:-20, from=0, to=.99, colors=ccol, ylim = c(.01, 1e7),
     log = "y", yaxt = "n")
if(require(sfsmisc)) eaxis(2) else axis(2)
}
\keyword{arith}
