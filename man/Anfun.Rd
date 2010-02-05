\name{Anfun}
\alias{Anfun}

\title{Nonparametric rank-based estimators of the Pickands dependence
  function in the bivariate case}
\description{
  Nonparametric rank-based estimators of the Pickands dependence
  function studied in Genest and Segers (2009). 
}
\usage{
Anfun(x, w, estimator = "CFG", corrected = TRUE)
}

\arguments{
  \item{x}{ a two-column data matrix that will be transformed to pseudo-observations. }
  \item{w}{ a vector of points in [0,1] where to evaluate the estimated
    Pickands dependence function.}
  \item{estimator}{ specifies which nonparametric rank-based estimator
    of the unknown Pickands dependence function to use; can be either
    \code{"CFG"} (Caperaa-Fougeres-Genest) or  \code{"Pickands"}.}
  \item{corrected}{ TRUE means that the estimators will be corrected to
    ensure that their value at 0 and 1 is 1.}
}
\details{
  More details can be found in Genest and Segers (2009).
}
\value{
  Returns a vector containing the values of the estimated
  Pickands dependence function at the points in \code{w}.
}
\references{

  C. Genest and J. Segers (2009). Rank-based inference for bivariate
  extreme-value copulas. \emph{Annals of Statistics}, 37, pages 2990-3022.
}
\seealso{ \code{\link{evTestC}}, \code{\link{evTestK}},
  \code{\link{evCopula}}, \code{\link{gofEVCopula}}. }
\examples{
## True Pickands dependence functions
curve(Afun(gumbelCopula(4), x), 0, 1)
curve(Afun(gumbelCopula(2), x), add=TRUE)
curve(Afun(gumbelCopula(1.33), x), add=TRUE)

## CFG estimator
curve(Anfun(rcopula(gumbelCopula(4), 1000), x), lty=2, add=TRUE)
curve(Anfun(rcopula(gumbelCopula(2), 1000), x),  lty=2, add=TRUE)
curve(Anfun(rcopula(gumbelCopula(1.33), 1000), x),  lty=2, add=TRUE)

## Pickands estimator
curve(Anfun(rcopula(gumbelCopula(4), 1000), x, estimator="Pickands"),
      lty=3, add=TRUE)
curve(Anfun(rcopula(gumbelCopula(2), 1000), x, estimator="Pickands"),
      lty=3, add=TRUE)
curve(Anfun(rcopula(gumbelCopula(1.33), 1000), x, estimator="Pickands"),
      lty=3, add=TRUE)
}

\keyword{htest}
\keyword{multivariate}
