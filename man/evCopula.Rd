\name{evCopula}
\alias{evCopula}
\alias{galambosCopula}
\alias{huslerReissCopula}
\alias{tawnCopula}
\alias{tevCopula}

\title{ Construction of extreme-value copula class objects }
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
    \code{df} will be estimated if the object is passed as argument to \code{fitCopula}.}

  \item{\dots}{ currently nothing. }
}
\value{
  An object of class \code{"gumbelCopula"}, \code{"galambosCopula"},
  \code{"huslerReissCopula"}, \code{"tawnCopula"}, or \code{"tevCopula"}.
}
\note{
  The Gumbel copula is both an Archimedean and an extreme-value copula.
}

\seealso{ \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{gofEVCopula}}, \code{\link{Anfun}}. }

\examples{
## Notice that, for a given degree of dependence,
## these copulas are strikingly similar.

tau <- 0.33

gumbel.cop <- evCopula("gumbel", param=2)
gumbel.cop@parameters <- calibKendallsTau(gumbel.cop, tau)

galambos.cop <- galambosCopula(2)
galambos.cop@parameters <- calibKendallsTau(galambos.cop, tau)

huslerReiss.cop <- huslerReissCopula(2)
huslerReiss.cop@parameters <- calibKendallsTau(huslerReiss.cop, tau)

tawn.cop <- tawnCopula(0.3)
tawn.cop@parameters <- calibKendallsTau(tawn.cop, tau)

tev.cop <- tevCopula(0.7)
tev.cop@parameters[1] <- calibKendallsTau(tev.cop, tau)

curve(Afun(gumbel.cop, x), 0, 1)
curve(Afun(galambos.cop, x), lty=2, add=TRUE)
curve(Afun(huslerReiss.cop, x), lty=3, add=TRUE)
curve(Afun(tawn.cop, x), lty=4, add=TRUE)
curve(Afun(tev.cop, x), lty=5, add=TRUE)
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
