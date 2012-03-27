\name{Bernoulli}
\alias{Bernoulli}
\alias{Bernoulli.all}
\title{Compute Bernoulli Numbers}
\description{
  Compute the \eqn{n}th Bernoulli number, or
  generate all Bernoulli numbers up to the \eqn{n}th,
  using diverse \code{method}s, that is, algorithms.

  \bold{NOTE} the current default methods will be changed -- to get
  better accuracy !
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
    applied. \code{"A-T"}, which stands for the Akiyama-Tanigawa algorithm
    may be nice and simple, but has bad numerical properties.
  \code{"sumRamanujan"} is somewhat more efficient but not yet
  implemented.}
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

system.time(B50  <- Bernoulli.all(50))# still "no time"
system.time(B100 <- Bernoulli.all(100))# still less than a milli second

## Using Bernoulli() is not much slower, but hopefully *more* accurate!
## Check first - TODO
system.time(B.1c <- Bernoulli(100))
## reset the cache:
assign("Bern.tab", list(), envir = copula:::.nacopEnv)

## BUT -- the algorithm is *really* not accurate enough ...
## ---> try to work with higher precision

## NB: The following does not print *unless* you evaluate it *outside*
##     the if(..) clause
if(require("Rmpfr")) { ## note that it has its own Bernoulli() !
 system.time(B100.250 <- as.numeric(Bernoulli.all(100, prec = 250)))
 ## 0.75 sec [Core i5 (2010)]
 m <- cbind(Bn = B100.250, "log10(rel.Err)" =
	    round(log1p( - B100/B100.250)/log(10), 2))
 rownames(m) <- paste("n=",0:100, sep="")
 m[1:5,]
 m[2*(1:15) -1,] ## for n=10: still 8 correct digits

 system.time(B100.1k <- as.numeric(Bernoulli.all(100, prec = 1024)))
 ## The first 34 are "the same", but after [41],
 ## even 250 precBits were *not* sufficient:
 round(log10(abs(1 - B100.250/B100.1k))[seq(1,99,by=2)], 2)

 ## some accuracy investigation:
 nn <- 12:80; nn <- nn[nn \%\% 2 == 0]; nn
 B.asy  <- sapply(nn, copula::Bernoulli, method="asymp")
 B.sumB <- sapply(nn, copula::Bernoulli, method="sumBin")
 B.prec <- Rmpfr::Bernoulli(nn, precBits = 512)
 relErr <- as.numeric(1 - B.asy  / B.prec)
 relE2 <-  as.numeric(1 - B.sumB / B.prec)
 matplot(nn, abs(cbind(relErr, relE2)),
	 ylim = c(1e-15, 1e-4), log="y", type="b")
 ##--> an optimal "hybrid" method will use "asymp" from about n ~= 20

}## end if(require("Rmpfr"))
}
\keyword{arithmetic}
