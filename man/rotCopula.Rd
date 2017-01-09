\name{rotCopula}
\title{Construction and Class of Rotated aka Reflected Copulas}
\concept{survival copula}
\concept{rotated copula}
\concept{reflected copula}
\alias{rotCopula}
\alias{rotCopula-class}
%
\alias{dim,rotCopula-method}
\description{
  Constructs a \dQuote{reflected} or \dQuote{rotated} copula from an
  initial copula and a vector of logicals indicating which dimension to
  \dQuote{flip}.
}
\usage{
rotCopula(copula, flip = TRUE)
}
\arguments{
  \item{copula}{an object of class \code{"\linkS4class{Copula}"}.}
  \item{flip}{\code{\link{logical}} vector (of length 1 or
    \code{dim(copula)}) indicating which dimension should be
    \dQuote{flipped}; by default, all the components are flipped,
    implying that the resulting copula is the \emph{\dQuote{survival copula}}.}
}
\value{
  A \dQuote{rotated} or \dQuote{reflected} copula object of class \code{"rotCopula"}.
}
\note{When there are an even number of flips, the resulting copula
  can be seen as a \emph{rotated} version of \code{copula}.  In the
  other cases, e.g., \code{flip = c(FALSE,TRUE)} in 2d, it is rather a
  a \emph{reflected} or \dQuote{mirrored} copula.
}
\section{Slots}{
  of a \code{"rotCopula"} object
  \describe{
    \item{\code{copula}:}{Object of class \code{"\linkS4class{copula}"}.}
    \item{\code{flip}:}{\code{\link{logical}} vector of length \eqn{d}
      (the copula dimension) specifying which margins are flipped;
      corresponds to the \code{flip} argument of \code{rotCopula()}. }
    \item{\code{dimension}:}{the copula dimension \eqn{d}, an \code{\link{integer}}.}
    \item{\code{parameters}:}{\code{\link{numeric}} vector specifying
      the parameters.}
    \item{\code{param.lowbnd}, and \code{param.upbnd}:}{\code{numeric}
      vector of the same length as \code{parameters}, specifying
      (component wise) bounds for each of the parameters.}
    \item{\code{param.names}:}{\code{\link{character}} vector (of same
      length as \code{parameters}) with parameter names.}
    \item{\code{fullname}:}{\bold{deprecated}; a \code{character} string describing the
      rotated copula.}
  }
}
%\details{}
%\author{Ivan Kojadinovic}
%\references{}
\seealso{
  \code{\link{fitCopula}} for fitting such copulas to data,
  \code{\link{gofCopula}} for goodness-of-fit tests for such copulas.
}
\examples{
## Two-dimensional examples: First the Clayton(3) survival copula:
survC <- rotCopula(claytonCopula(3)) # default: flip = 'all TRUE'
contourplot2(survC, dCopula)

## Now, a reflected Clayton copula:
r10C  <- rotCopula(claytonCopula(3), flip = c(TRUE, FALSE))

contourplot2(r10C, dCopula, nlevels = 20, main = "dCopula(<rotCopula>)")
contourplot2(r10C, pCopula, nlevels = 20, main = "pCopula(<rotCopula>)")
rho(r10C)
tau(r10C) # -0.6

n <- 1000
u <- rCopula(n, r10C)
rho.n <- cor(u[,1], u[,2], method = "spearman")
tau.n <- cor(u[,1], u[,2], method = "kendall")

## "Calibration"
rc. <- rotCopula(claytonCopula(), flip = c(TRUE, FALSE))
iRho(rc., rho.n)
iTau(rc., tau.n)

## Fitting
fitCopula(rc., data = pobs(u), method = "irho")
fitCopula(rc., data = pobs(u), method = "itau")
fitCopula(rc., data = pobs(u), method = "mpl")

## Goodness-of-fit testing -- the first, parametric bootstrap, is *really* slow
\dontrun{gofCopula(rc., x = u)}
\donttest{gofCopula(rc., x = u, simulation = "mult")}

## A four-dimensional example: a rotated Frank copula
rf <- rotCopula(frankCopula(10, dim = 4),
                flip = c(TRUE, FALSE, TRUE, FALSE))

n <- 1000
u <- rCopula(n, rf)
splom2(u)

pCopula(c(0.6,0.7,0.6,0.8), rf)
C.n(cbind(0.6,0.7,0.6,0.8),  u)

## Fitting
(rf. <- rotCopula(frankCopula(dim=4),
                  flip = c(TRUE, FALSE, TRUE, FALSE)))
## fitCopula(rf., data = pobs(u), method = "irho")
## FIXME above: not related to rotCopula but frankCopula
fitCopula(rf., data = pobs(u), method = "itau")
fitCopula(rf., data = pobs(u), method = "mpl")

## Goodness-of-fit testing (first ("PB") is really slow, see above):
\dontrun{gofCopula(rf., x = u)}
\donttest{gofCopula(rf., x = u, simulation = "mult") # takes 3.7 sec [lynne, 2015]}
}
\keyword{distribution}
\keyword{multivariate}
