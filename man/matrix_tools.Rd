\name{matrix_tools}
\title{Tools to Work with Matrices}
\alias{p2P}
\alias{P2p}
\alias{getSigma}
\alias{extremePairs}
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

  \code{extremePairs()} identifies pairs with the largest (or smallest
  or both) entries in a symmetric matrix.
}
\usage{
p2P(param, d = floor(1 + sqrt(2*length(param))))
P2p(P)
getSigma(copula)
extremePairs(x, n = 6, method = c("largest", "smallest", "both"),
             use.names = FALSE)
}
\arguments{
  \item{param}{a parameter vector.}
  \item{d}{dimension of the resulting \code{\link{matrix}}.  The default
  is correct under the assumption (of \code{p2P()} in general!)  that
  \code{param} is the lower-triangular part of a correlation matrix
  \code{P} and hence corresponds to \code{\link{ellipCopula}(.., dispstr
  = "un")}.}
  \item{P}{a \code{\link{matrix}} which should be converted to a vector.}

  \item{copula}{an \bold{elliptical} copula, i.e., an object (extending)
    class \code{\linkS4class{ellipCopula}}; typically resulting from
    \code{\link{tCopula}()} or \code{\link{normalCopula}()}.}
  \item{x}{a symmetric \code{\link{matrix}}.}
  \item{n}{the number of pairs with smallest (or largest) values to be
    displayed.}
  \item{method}{a \code{\link{character}} string indicating the method
    to be used (with \code{"largest"} to comute the \code{n} pairs with
    largest entries in \code{x} (sorted in decreasing order); with
    \code{"smallest"} to compute the \code{n} pairs with smallest entries
    in \code{x} (sorted in increasing order); and with \code{"both"} to
    comute the 2\code{n} pairs with \code{n} largest entries and
    \code{n} smallest entries (sorted in decreasing order)).}
  \item{use.names}{A \code{\link{logical}} indicating whether
    \code{colnames(x)} are used as labels (if \code{!is.null(colnames(x))}).}
}
\value{
  \describe{
    \item{\code{p2P}:}{a symmetric \code{\link{matrix}} with ones on the
      diagonal and the values of \code{param} filled column-wise below the
      diagonal (which corresponds to row-wise filling above the diagonal).}
    \item{\code{P2p}:}{\code{\link{vector}} of column-wise below-diagonal
      entries of \code{P} (equal to the row-wise above-diagonal entries in
      case of a symmetric matrix).}
    \item{\code{getSigma}:}{\code{\link{matrix}} as from \code{p2P()}
      for all cases of elliptical copulas.}
    \item{\code{extremePairs}:}{a \code{\link{data.frame}}
      consisting of three columns (\code{row} (index or name), \code{col} (index or name),
      \code{value}).}
  }
}
\details{
  These auxiliary functions are often used when working with
  elliptical copulas.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{tCopula}}, \code{\link{normalCopula}}.
}
%\author{Marius Hofert and Martin M}
\examples{
## display the definitions
p2P
P2p
extremePairs

param <- (2:7)/10
tC <- tCopula(param, dim = 4, dispstr = "un", df = 3)
## consistency of the three functions :
P <- p2P(param) # (using the default 'd')
stopifnot(identical(param, P2p(P)),
	  identical(P, getSigma(tC)))

## Toeplitz case:
(tCt <- tCopula((2:6)/10, dim = 6, disp = "toep"))
(rhoP <- tCt@getRho(tCt))
stopifnot(identical(getSigma (tCt),
		    toeplitz (c(1, rhoP))))

## "AR1" case:
nC.7 <- normalCopula(0.8, dim = 7, dispstr = "ar1")
(Sar1.7 <- getSigma(nC.7))
0.8^(0:(7-1)) #  1  0.8  0.64  0.512 ..
stopifnot(all.equal(Sar1.7, toeplitz(0.8^(0:(7-1)))))
}
\keyword{array}
\keyword{manip}

