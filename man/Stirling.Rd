\name{Stirling}
\title{Eulerian and Stirling numbers of first and second kind}
\alias{Eulerian}
\alias{Stirling1}
\alias{Stirling2}
\alias{Eulerian.all}
\alias{Stirling1.all}
\alias{Stirling2.all}
\description{
  Compute Eulerian numbers and Stirling numbers of the first and second
  kind, possibly vectorized for all \eqn{k} \dQuote{at once}.
}
\usage{
Stirling1(n, k)
Stirling2(n, k, method = c("lookup.or.store", "direct"))
Eulerian (n, k, method = c("lookup.or.store", "direct"))

Stirling1.all(n)
Stirling2.all(n)
Eulerian.all (n)
}
\arguments{
  \item{n}{positive integer (\code{0} is allowed for \code{Eulerian()}).}
  \item{k}{integer in \code{0:n}.}
  \item{method}{for \code{Eulerian()} and \code{Stirling2()}, string
    specifying the method to be used.  \code{"direct"} uses the explicit
    formula (which may suffer from some cancelation for \dQuote{large}
    \code{n}).}
}
\details{
  Eulerian numbers:\cr
  \eqn{A(n,k) =} the number of permutations of 1,2,\dots,n with exactly \eqn{k}
  ascents (or exactly \eqn{k} descents).

  Stirling numbers of the first kind:\cr
  s(n,k) = (-1)^{n-k} times
  the number of permutations of 1,2,\dots,n with exactly k cycles.

  Stirling numbers of the second kind:\cr
  \eqn{S^{(k)}_n}{S(n,k)} is the number of ways of partitioning a set
    of \eqn{n} elements into \eqn{k} non-empty subsets.
}
\value{
  \eqn{A(n,k)}, \eqn{s(n,k)} or \eqn{S(n,k) = S^{(k)}_n}{S(n,k)}, respectively.

  \code{Eulerian.all(n)} is the same as \code{sapply(0:(n-1), Eulerian, n=n)}
  (for \eqn{n > 0}), \cr
  \code{Stirling1.all(n)} is the same as \code{sapply(1:n, Stirling1, n=n)},
  and\cr
  \code{Stirling2.all(n)} is the same as \code{sapply(1:n, Stirling2, n=n)},
  but more efficient.
}
\note{
  For typical double precision arithmetic,\cr
  \code{Eulerian*(n, *)} overflow (to \code{Inf}) for \eqn{n \ge 172},\cr
  \code{Stirling1*(n, *)} overflow (to \eqn{\pm}{+/-}\code{Inf}) for
  \eqn{n \ge 171}, and\cr
  \code{Stirling2*(n, *)} overflow (to \code{Inf}) for \eqn{n \ge 220}.
%% and it would be possible to implement Stirling*(n, log=TRUE) without
%% overflow, but a bit tedious, as one should "log-scale arithmetic only
%% for large n, and hence as to *switch* for n >= n.max

}
\author{Martin Maechler ("direct": May 1992)}
\references{
  \bold{Eulerians:}

  NIST Digital Library of Mathematical Functions,
  26.14: \url{http://dlmf.nist.gov/26.14}

  \bold{Stirling numbers:}

  Abramowitz and Stegun
  24,1,4 (p. 824-5 ; Table 24.4, p.835);
  Closed Form : p.824 "C."

  NIST Digital Library of Mathematical Functions,
  26.8: \url{http://dlmf.nist.gov/26.8}
}
% \seealso{
% }
\examples{
Stirling1(7,2)
Stirling2(7,3)

Stirling1.all(9)
Stirling2.all(9)
}
\keyword{arithmetic}
