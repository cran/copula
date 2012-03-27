\name{safeUroot}
\alias{safeUroot}
\title{One Dimensional Root (Zero) Finding - Extra "Safety" for Convenience}
\usage{
safeUroot(f, interval, \dots,
       lower = min(interval), upper = max(interval),
       f.lower = f(lower, \dots), f.upper = f(upper, \dots),
       Sig = NULL, check.conv = FALSE,
       tol = .Machine$double.eps^0.25, maxiter = 1000, trace = 0)
}
\description{
  \code{safeUroot()} as a \dQuote{safe} version of
  \code{\link{uniroot}()} searches for a root (that is, zero) of the
  function \code{f} with respect to its first argument.

  \dQuote{Safe} means searching for the correct \code{interval
    = c(lower,upper)} if \code{sign(f(x))} does not satisfy the
  requirements at the interval end points; see the \sQuote{Details}
  section.
}
\arguments{
  \item{f}{function}
  \item{interval}{interval}
  \item{\dots}{additional named or unnamed arguments to be passed
    to \code{f}}
  \item{lower, upper}{lower and upper endpoint of search interval}
  \item{f.lower, f.upper}{function value at \code{lower} or \code{upper}
    endpoint, respectively.}
  \item{Sig}{\emph{desired} sign of \code{f(upper)}, or \code{\link{NULL}}.}
  \item{check.conv}{logical indicating whether a convergence warning of the
    underlying \code{\link{uniroot}} should be caught as an error.}
  \item{tol}{the desired accuracy, that is, convergence tolerance.}
  \item{maxiter}{maximal number of iterations}
  \item{trace}{number determining tracing}
}
\details{
  If it is \emph{known how} \eqn{f} changes sign at the root
  \eqn{x_0}{x0}, that is, if the function is increasing or decreasing there,
  \code{Sig} can be specified, typically as \eqn{S := \pm 1}{S:= +/- 1},
  to require \eqn{S = \mathrm{sign}(f(x_0 + \epsilon))}{S = sign(f(x0 + eps))}
  at the solution.  In that case, the search interval \eqn{[l,u]} must
  be such that \eqn{S * f(l) <= 0} and \eqn{S * f(u) >= 0}.

  Otherwise, by default, when \code{Sig=NULL}, the search interval
  \eqn{[l,u]} must satisfy \eqn{f(l) * f(u) <= 0}.

  In both cases, when the requirement is not satisfied, \code{safeUroot()}
  tries to enlarge the interval until the requirement \emph{is} satisfied.

  %% in the Sig case, typically only *one* side of the interval is
  %% enlarged; in the other (usual) case, it is enlarged on both sides
  %% at the same time.
}
\value{
  A list with four components, \code{root}, \code{f.root}, \code{iter}
  and \code{estim.prec}; see \code{\link{uniroot}}.
}
\author{Martin Maechler (from Martin's \R package \pkg{nor1mix}).
}
\seealso{
  \code{\link{uniroot}}.
}
\examples{
f1 <- function(x) (121 - x^2)/(x^2+1)
f2 <- function(x) exp(-x)*(x - 12)

try(uniroot(f1, c(0,10)))
try(uniroot(f2, c(0,2)))
##--> error: f() .. end points not of opposite sign

## where as safeUroot() simply first enlarges the search interval:
safeUroot(f1, c(0,10),trace=1)
safeUroot(f2, c(0,2), trace=2)

## no way to find a zero of a positive function:
try( safeUroot(exp, c(0,2), trace=TRUE) )

## Convergence checking :
safeUroot(sinc, c(0,5), maxiter=4) #-> "just" a warning
try( # an error, now with  check.conv=TRUE
  safeUroot(sinc, c(0,5), maxiter=4, check.conv=TRUE) )
}
\keyword{optimize}
