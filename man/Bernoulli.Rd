\name{Bernoulli}
\alias{Bernoulli}
\alias{Bernoulli.all}
\title{Compute Bernoulli Numbers}
\description{
  Compute the \eqn{n}th Bernoulli number, or
  generate all Bernoulli numbers up to the \eqn{n}th,
  using diverse \code{method}s, that is, algorithms.

  \bold{NOTE} the current default methods will be changed -- to get
  better accuracy!
%%--> ../R/special-func.R
}
\usage{
Bernoulli    (n, method = c("sumBin", "sumRamanujan", "asymptotic"),
              verbose = FALSE)
Bernoulli.all(n, method = c("A-T", "sumBin", "sumRamanujan", "asymptotic"),
              precBits = NULL, verbose = getOption("verbose"))
}
\arguments{
  \item{n}{positive integer, indicating the index of the largest (and
    last) of the Bernoulli numbers needed.}
  \item{method}{character string, specifying which method should be
    applied. The default for \code{Bernoulli.all()}, \code{"A-T"}
    stands for the Akiyama-Tanigawa algorithm
    which is nice and simple but has bad numerical properties.  It can
    however work with high precision "mpfr"-numbers, see \code{precBits}.
    \code{"sumRamanujan"} is somewhat more efficient but not yet implemented.}
  \item{precBits}{currently only for \code{method = "A-T"} --
    \code{NULL} or a positive integer indicating the precision of the
    initial numbrs in bits, using \pkg{"Rmpfr"}'s package multiprecision
    arithmetic.}
  \item{verbose}{(for \code{"A-T"}:) logical indicating if the
    intermediate results of the algorithm should be printed.}
}
\value{
  \describe{
    \item{\code{Bernoulli()}:}{a number}
    \item{\code{Bernoulli.all()}:}{a numeric vector of length n,
      containing B(n)}
  }
}
\author{Martin Maechler}% started: 25 Jun 2011 (night train Vienna - Zurich).
% \details{
% The usual algorithms use O(n^2) math operations and in addition lose
% precision rapidly.
% }
\references{
  Kaneko, Masanobu (2000)
  The Akiyama-Tanigawa algorithm for Bernoulli numbers;
  Journal of Integer Sequences \bold{3}, article 00.2.9
}
\seealso{
  \code{\link{Eulerian}}, \code{\link{Stirling1}}, etc.
}
\examples{
## The example for the paper
MASS::fractions(Bernoulli.all(8, verbose=TRUE))

B10 <- Bernoulli.all(10)
MASS::fractions(B10)

system.time(B50  <- Bernoulli.all(50))#  {does not cache} -- still "no time"
system.time(B100 <- Bernoulli.all(100))# still less than a milli second

## Using Bernoulli() is not much slower, but hopefully *more* accurate!
## Check first - TODO
system.time(B.1c <- Bernoulli(100))# caches ..
system.time(B1c. <- Bernoulli(100))# ==> now much faster
stopifnot(identical(B.1c, B1c.))

if(FALSE)## reset the cache:
assign("Bern.tab", list(), envir = copula:::.nacopEnv)

## More experiments in the source of the copula package ../tests/Stirling-etc.R
}
\keyword{arith}
