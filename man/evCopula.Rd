\name{evCopula}
\alias{evCopula}
\alias{galambosCopula}
\alias{huslerReissCopula}
\alias{tawnCopula}
\alias{tevCopula}
\title{Construction of Extreme-Value Copula Objects}
\description{
  Constructs an extreme-value copula class object with its corresponding
  parameter.
}
\usage{
evCopula(family, param, dim = 2, ...)
galambosCopula(param)
huslerReissCopula(param)
tawnCopula(param)
tevCopula(param, df = 4, df.fixed = FALSE)
%gumbelCopula(param) ==> ./archmCopula.Rd
}
\arguments{
  \item{family}{a character string specifying the family of an
    extreme-value copula.  Currently implemented are %% <--> ../R/evCopula.R
    \code{"galambos", "gumbel", "huslerReiss", "tawn", "tev"}.}
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{the dimension of the copula.  Currently, only \code{"gumbel"}
    allows \code{dim > 2}.}
  \item{df}{a number specifying the degrees of freedom of the t
    extreme-value copula, \code{tevCopula()}.}
  \item{df.fixed}{logical; if true, the degrees of freedom will never be
    considered as a parameter to be estimated; \code{FALSE} means that
    \code{df} will be estimated if the object is passed as argument to
    \code{\link{fitCopula}}.}

  \item{\dots}{ currently nothing. }
}
\value{
  An object of class \code{"\linkS4class{gumbelCopula}"},
  \code{"\linkS4class{galambosCopula}"},
  \code{"\linkS4class{huslerReissCopula}"},
  \code{"\linkS4class{tawnCopula}"}, or \code{"\linkS4class{tevCopula}"}.
}
\note{
  The Gumbel copula is both an Archimedean and an extreme-value copula,
  with principal documentation on \code{\link{gumbelCopula}} (or
  \code{\link{archmCopula}}).
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{ \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}. }
\examples{
## Gumbel is both
stopifnot(identical(   evCopula("gumbel"), gumbelCopula()),
          identical(archmCopula("gumbel"), gumbelCopula()))

## For a given degree of dependence these copulas are strikingly similar :

tau <- 1/3

gumbel.cop      <- gumbelCopula     (iTau(gumbelCopula(),      tau))
galambos.cop    <- galambosCopula   (iTau(galambosCopula(),    tau))
huslerReiss.cop <- huslerReissCopula(iTau(huslerReissCopula(), tau))
tawn.cop        <- tawnCopula       (iTau(tawnCopula(),        tau))
tev.cop         <- tevCopula        (iTau(tevCopula(),         tau))

curve(A(gumbel.cop, x), 0, 1, ylab = "A(<cop>( iTau(<cop>(), tau)), x)",
      main = paste("A(x) for five Extreme Value cop. w/  tau =", format(tau)))
curve(A(galambos.cop, x), lty=2, add=TRUE)
curve(A(huslerReiss.cop, x), lty=3, add=TRUE)
curve(A(tawn.cop, x), lty=4, add=TRUE)
curve(A(tev.cop, x), lty=5, col=2, add=TRUE)# very close to Gumbel

## And look at the differences
curve(A(gumbel.cop, x) - A(tawn.cop, x), ylim = c(-1,1)*0.005,
      ylab = '', main = "A(<Gumbel>, x) - A(<EV-Cop.>, x)")
abline(h=0, lty=2)
curve(A(gumbel.cop, x) - A(galambos.cop, x), add=TRUE, col=2)
curve(A(gumbel.cop, x) - A(huslerReiss.cop, x), add=TRUE, col=3)
curve(A(gumbel.cop, x) - A(tev.cop, x), add=TRUE, col=4, lwd=2)


## the t-EV-copula has always positive tau :
curve(vapply(x, function(x) tau(tevCopula(x)), 0.), -1, 1,
      n=257, ylim=0:1, xlab=quote(rho),ylab=quote(tau),
      main= quote(tau( tevCopula(rho) )), col = 2, lwd = 2)
rect(-1,0,1,1, lty = 2, border = adjustcolor("black", 0.5))
}
\keyword{distribution}
\keyword{multivariate}
