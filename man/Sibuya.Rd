\name{Sibuya}
\alias{rSibuya}
\alias{pSibuya}
\alias{dSibuya}
\alias{rSibuyaR}
\alias{dsumSibuya}
\title{Sibuya distribution - sampling and probabilities}
\description{
  The Sibuya distribution \eqn{\mathrm{Sib}(\alpha)}{Sib(alpha)} can be
  defined by its Laplace transform
  \deqn{1-(1-\exp(-t))^\alpha,\ t\in[0,\infty),}{1-(1-exp(-t))^alpha, t
 in [0,Inf),}
  its distribution function% "R man - FIXME: allow to use \mathbb{N} instead:
  \deqn{F(k)=1-(-1)^k{\alpha-1\choose k}=1-\frac{1}{kB(k,1-\alpha)},\
    k\in\mathbf{N}}{%
    F(k)=1-(-1)^k * choose(alpha-1,k)=1-1/(k * beta(k, 1-alpha)), k in N,}
  (where \eqn{B}{beta} denotes the beta function) or its probability
  mass function \deqn{p_k={\alpha\choose k}(-1)^{k-1},\ k\in\mathbf{N},
  }{p_k = choose(alpha,k)*(-1)^(k-1), k in N,}
  where \eqn{\alpha\in(0,1]}{alpha in (0,1]}.

  \code{pSibuya} evaluates the distribution function.

  \code{dSibuya} evaluates the probability mass function.

  \code{rSibuya} generates random variates from
  \eqn{\mathrm{Sib}(\alpha)}{Sib(alpha)} with
  the algorithm described in Hofert (2011), Proposition 3.2.

  \code{dsumSibuya} gives the probability mass function of the
  \eqn{n}-fold convolution of Sibuya variables, that is, the sum of \eqn{n}
  independent Sibuya random variables,
  \eqn{S = \sum_{i=1}^n X_i}{S = sum(i=1..n) X[i]}, where
  \eqn{X_i \sim \mathrm{Sib}(\alpha)}{X[i] ~ Sib(alpha)}.

  This probability mass function can be shown (see Hofert
  (2010, pp. 99)) to be
  \deqn{\sum_{j=1}^n{n\choose j}{j\alpha\choose k} (-1)^{k-j},\
    k\in\{n,n+1,\dots\}.}{%
    sum(j=1..n) choose(n,j) * choose(j*alpha,k) * (-1)^(k-j), k >= n.}
}
\usage{
rSibuya(n, alpha)
dSibuya(x, alpha, log=FALSE)
pSibuya(x, alpha, lower.tail=TRUE, log.p=FALSE)

dsumSibuya(x, n, alpha,
           method=c("log", "direct", "diff", "exp.log",
                    "Rmpfr", "Rmpfr0", "RmpfrM", "Rmpfr0M"),
           mpfr = list(minPrec = 21, fac = 1.25, verbose=TRUE),
           log=FALSE)
}
\arguments{
  \item{n}{
    for \code{rSibuya}: sample size, that is, length of the resulting
    vector of random variates.
    \cr
    for \code{dsumSibuya}: the number \eqn{n} of summands.
  }
  \item{alpha}{parameter in \eqn{(0,1]}.}
  \item{x}{vector of \code{\link{integer}} values (\dQuote{quantiles})
\eqn{x} at which to compute the probability mass or cumulative probability.}
  \item{log, log.p}{\code{\link{logical}}; if TRUE, probabilities p are
given as log(p).}
  \item{lower.tail}{\code{\link{logical}}; if TRUE (the default), probabilities
      are \eqn{P(X \le x)}{P(X <= x)}, otherwise, \eqn{P(X > x)}.}
  \item{method}{character string specifying which computational method is to be
    applied.  Implemented are:
    \describe{
      \item{\code{"log"}}{evaluates the logarithm of the sum
	\deqn{\sum_{j=1}^n {n\choose j}{j\alpha\choose x} (-1)^{x-j}}{%
	  sum(j=1..n) choose(n,j) * choose(j*alpha,x) * (-1)^(x-j)}
	in a numerically stable way;}
      \item{\code{"direct"}}{directly evaluates the sum;}
      \item{\code{"Rmpfr*"}}{are as \code{method="direct"} but use
        high-precision arithmetic; \code{"Rmpfr"} and \code{"Rmpfr0"} return
        \code{\link{double}}s whereas \code{"RmpfrM"} and \code{"Rmpfr0M"} give
        \code{\link[Rmpfr]{mpfr}} high-precision numbers.
	Whereas \code{"Rmpfr"} and \code{"RmpfrM"} each adapt to high
        enough precision, the \code{"Rmpfr0*"} ones do not adapt.\cr
%% FIXME --- when do we *compute* the precision (without adaption, but
%% by pre() formula? --- exactly when 'alpha' is not an mpfr number
	For all \code{"Rmpfr*"} methods, \code{alpha} can be set to a
        \code{\link[Rmpfr:mpfr-class]{mpfr}} number of specified
        precision and this will determine the precision of all parts of
    the internal computations.
      }
      \item{\code{"diff"}}{interprets the sum as a forward difference
	and computes it via \code{diff};}
      \item{\code{"exp.log"}}{is as \code{method="log"} but without
	numerically stable evaluation (not recommended, use with care).}
    }
  }
  \item{mpfr}{for \code{method = "Rmpfr"} or \code{"RmpfrM"} only: a
    list of
    \cr \code{minPrec}: minimal (estimated) precision in bits,
    \cr \code{fac}: factor with which current precision is multiplied if
    it is not sufficient.
    \cr \code{verbose}: determining if and how much is printed.
  }
}
\details{
  The Sibuya distribution has \bold{no} finite moments, that is, specifically
  infinite mean and variance.

  For documentation and didactical purposes, \code{rSibuyaR} is a pure-\R
  implementation of \code{rSibuya}, of course slower than \code{rSibuya}
  as the latter is implemented in C.

  Note that the sum to evaluate for \code{dsumSibuya} is numerically
  highly challenging, even already for small
  \eqn{\alpha}{alpha} values (for example, \eqn{n \ge 10}{n >= 10}),
  and therefore should be used with care.  It may require high-precision
  arithmetic which can be accessed with \code{method="Rmpfr"} (and the
  \pkg{Rmpfr} package).
}
\value{
  \describe{
    \item{rSibuya:}{A vector of positive \code{\link{integer}}s of
      length \code{n} containing the generated random variates.}
    \item{dSibuya, pSibuya:}{a vector of
      probabilities of the same length as \code{x}.}
    \item{dsumSibuya:}{a vector of probabilities, positive if and only if
      \code{x >= n} and of the same length as \code{x} (or \code{n} if
      that is longer).}
  }
}
\author{Marius Hofert, Martin Maechler}
\references{
  Hofert, M. (2010).
  \emph{Sampling Nested Archimedean Copulas with Applications to CDO Pricing}.
  S\enc{ü}{ue}dwestdeutscher Verlag fuer Hochschulschriften AG & Co. KG.

  Hofert, M. (2011).
  Efficiently sampling nested Archimedean copulas.
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.
}
\seealso{
  \code{\link{rFJoe}} and \code{\link{rF01Joe}} (where \code{rSibuya} is
  applied).
}
\examples{
## Sample n random variates from a Sibuya(alpha) distribution and plot a
## histogram
n <- 1000
alpha <- .4
X <- rSibuya(n, alpha)
hist(log(X), prob=TRUE); lines(density(log(X)), col=2, lwd=2)
}
\keyword{distribution}
