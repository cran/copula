\name{evCopula}
\alias{evCopula}
\alias{galambosCopula}
\alias{huslerReissCopula}
\alias{tawnCopula}
\alias{tevCopula}
\title{Construction of extreme-value copula class objects}
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
%gumbelCopula(param)
}
\arguments{
  \item{family}{ a character string specifying the family of an
    extreme-value copula. %Implemented families are "gumbel", "galambos" and
    %"gumbel".
  }
  \item{param}{ a numeric vector specifying the parameter values. }
  \item{dim}{ the dimension of the copula. }
  \item{df}{ a numerical value specifying the number of degrees of
    freedom the t extreme-value copula.}
  \item{df.fixed}{ TRUE means that the degrees of freedom will never be
    considered as a parameter to be estimated; FALSE means that
    \code{df} will be estimated if the object is passed as argument to
    \code{fitCopula}.}

  \item{\dots}{ currently nothing. }
}
\value{
  An object of class \code{"\linkS4class{gumbelCopula}"},
  \code{"\linkS4class{galambosCopula}"},
  \code{"\linkS4class{huslerReissCopula}"},
  \code{"\linkS4class{tawnCopula}"}, or \code{"\linkS4class{tevCopula}"}.
}
\note{
  The Gumbel copula is both an Archimedean and an extreme-value copula.
}
\seealso{ \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}. }
\examples{
## Notice that, for a given degree of dependence,
## these copulas are strikingly similar.

tau <- 0.33

gumbel.cop <- evCopula("gumbel")
stopifnot(identical(gumbel.cop, gumbelCopula()))
gumbel.cop@parameters <- iTau(gumbel.cop, tau)

galambos.cop <- galambosCopula()
galambos.cop@parameters <- iTau(galambos.cop, tau)

huslerReiss.cop <- huslerReissCopula()
huslerReiss.cop@parameters <- iTau(huslerReiss.cop, tau)

tawn.cop <- tawnCopula()
tawn.cop@parameters <- iTau(tawn.cop, tau)

tev.cop <- tevCopula()
tev.cop@parameters[1] <- iTau(tev.cop, tau)

curve(A(gumbel.cop, x), 0, 1,
      main = "A(x) for five Extreme Value cop. w/  tau = 1/3")
curve(A(galambos.cop, x), lty=2, add=TRUE)
curve(A(huslerReiss.cop, x), lty=3, add=TRUE)
curve(A(tawn.cop, x), lty=4, add=TRUE)
curve(A(tev.cop, x), lty=5, add=TRUE)

## the t-EV-copula has always positive tau :
curve(vapply(x, function(x) tau(tevCopula(x)), 0.), -1, 1,
      col=2, n=257, ylim=0:1, ylab=bquote(tau), xlab=bquote(rho),
      main= expression(tau( tevCopula(rho) )))
}
\keyword{ distribution }
\keyword{ multivariate }
