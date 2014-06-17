\name{contour-methods}
\docType{methods}
\alias{contour-methods}
\alias{contour,indepCopula-method}
\alias{contour,copula-method}
\alias{contour,mvdc-method}
\title{Methods for Function 'contour' in Package 'copula'}
\description{
  Methods for function \code{\link{contour}} to draw contour lines aka a
  level plot for objectsf from package \pkg{copula}.
}
\usage{% want \usage{} since we have many "surprising arguments"
\S4method{contour}{copula}(x, fun,
      n = 51, delta = 0, box01 = TRUE, \dots)

\S4method{contour}{mvdc}(x, fun,
      xlim, ylim, nx = 51, ny = 51,
      xis = seq(xlim[1], xlim[2], length = nx),
      yis = seq(ylim[1], ylim[2], length = ny),
      box01 = TRUE, \dots)
}
\arguments{
  \item{x}{either a \code{"\linkS4class{copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{fun}{the \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}

  %% "copula" only :
  \item{n}{(for \code{"copula"}:) the number of points in \emph{both}
    directions to do the plotting.  The function \code{fun} will be
    evaluated on a grid of size \eqn{n \times n}{n x n}.}
  \item{delta}{a very small number in \eqn{[0, \frac 1 2)}{[0, 1/2)},
    defaulting to zero.  The x- and y-ranges to be used for plotting
    will be \code{[0+delta, 1-delta]}, i.e., \code{[0,1]} by default.}

  %% "mvdc" only :
  \item{xlim, ylim}{(\code{"mvdc"}:) the \code{\link{range}} of the x or
    y variable, respectively.}
  \item{nx,ny}{(\code{"mvdc"}:) the number of points in x- or
    y-direction, respectively.  The function \code{fun} will be
    evaluated on a grid of size \eqn{nx \times ny}{nx x ny}.}
  \item{xis, yis}{(\code{"mvdc"}:) instead of specifying \code{xlim,
      ylim} and \code{nx, ny}, the numeric vectors (of length \code{nx}
    and \code{ny}) may be specified directly.}

  \item{box01}{logical specifying if faint rectangle should be drawn
    at the \eqn{[0,1]^2} borders.}
  %% contour
  \item{\dots}{further arguments for (the default method of)
    \code{\link{contour}()}, e.g., \code{nlevels}, \code{levels}, etc.}
}
\section{Methods}{
  Contour lines are drawn for \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"} objects, see \code{x} in the
  \emph{Arguments} section.
}
\seealso{
  The \code{\link{persp-methods}} for \dQuote{perspective} aka
  \dQuote{3D} plots.
}
\examples{
contour(frankCopula(-0.8), dCopula)
contour(frankCopula(-0.8), dCopula, delta=1e-6)
contour(frankCopula(-1.2), pCopula)
contour(claytonCopula(2), pCopula)

## the Gumbel copula density is "extreme"
## --> use fine grid (and enough levels):
r <- contour(gumbelCopula(3), dCopula, n=200, nlevels=100)
range(r$z)# [0, 125.912]
## Now superimpose contours of three resolutions:
contour(r, levels=seq(1, max(r$z), by=2), lwd=1.5)
contour(r, levels = (1:13)/2, add=TRUE, col=adjustcolor(1,3/4), lty=2)
contour(r, levels = (1:13)/4, add=TRUE, col=adjustcolor(2,1/2),
        lty=3, lwd=3/4)

x <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd =1), list(mean = 1)))
contour(x, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
contour(x, pMvdc, xlim=c(-2, 2), ylim=c(-1, 3))
}
\keyword{methods}
\keyword{hplot}
