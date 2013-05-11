\name{p2P}
\alias{p2P}
\alias{P2p}
\alias{getSigma}
\title{Convert (Rho) Matrices to and From Parameter Vectors}
\description{
  \code{p2P()} creates a \code{\link{matrix}} from a given
  \code{\link{vector}} of parameters.  \code{P2p()} creates a
  numeric vector from a given \code{\link{matrix}},
  currently useful for elliptical copulas.

  \code{getSigma()} returns the \eqn{d \times d}{d * d} symmetric
  matrix \eqn{\Sigma} which is called \dQuote{Rho} as well, written
  (capital Greek \eqn{\rho}{rho} !) as \eqn{P} (and hence sometimes
  erronously pronounced "Pee").
  Note that \code{getSigma()} works for all elliptical copulas and uses
  \code{p2P()} for the \dQuote{unstuctured} case, \code{dispstr = "un"}.
}
\usage{
p2P(param, d)
P2p(P)
getSigma(copula)
}
\arguments{
  \item{param}{a parameter vector.}
  \item{d}{dimension of the resulting \code{\link{matrix}}.}
  \item{P}{\code{\link{matrix}} which should be converted to a vector.}
  \item{copula}{an \bold{elliptical} copula, i.e., an object (extending)
    class \code{\linkS4class{ellipCopula}}; typically resulting from
    \code{\link{tCopula}()} or \code{\link{normalCopula}()}.}
}
\value{
    \describe{
    \item{\code{p2P}: }{a symmetric \code{\link{matrix}} with ones on the
      diagonal and the values of \code{param} filled column-wise below the
      diagonal (which corresponds to row-wise filling above the diagonal).}
    \item{\code{P2p}: }{\code{\link{vector}} of column-wise below-diagonal
      entries of \code{P} (equal to the row-wise above-diagonal entries in
      case of a symmetric matrix).}
    \item{\code{getSigma}: }{\code{\link{matrix}} as from \code{p2P()}
      for all cases of elliptical copulas.}
  }
}
\details{
  These auxiliary functions are often used when working with
  elliptical copulas.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{tCopula}}, \code{\link{normalCopula}}.
}
%\author{Marius Hofert.}
\examples{
## display the definitions
p2P
P2p

param <- (2:7)/10
tC <- tCopula(param, dim = 4, dispstr="un", df = 3)
## consistency of the three functions :
P <- p2P(param, d=4)
stopifnot(identical(param, P2p(P)),
	  identical(P, getSigma(tC)))

## Toeplitz case:
(tCt <- tCopula((2:6)/10, dim = 6, disp = "toep"))
(rhoP <- tCt@getRho(tCt))
stopifnot(identical(getSigma (tCt),
		    toeplitz (c(1, rhoP))))
}
\keyword{array}
\keyword{manip}

