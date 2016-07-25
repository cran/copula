\name{contour-methods}
\title{Methods for Contour Plots in Package 'copula'}
\docType{methods}
\alias{contour-methods}
\alias{contour,indepCopula-method}
\alias{contour,Copula-method}
\alias{contour,mvdc-method}
\description{
  Methods for function \code{\link{contour}} to draw contour lines aka a
  level plot for objects from package \pkg{copula}.
}
\usage{% want \usage{} since we have many "surprising arguments"
\S4method{contour}{Copula}(x, FUN,
                   n.grid = 26, delta = 0,
                   xlab = expression(u[1]), ylab = expression(u[2]),
                   box01 = TRUE, \dots)
\S4method{contour}{mvdc}(x, FUN, xlim, ylim, n.grid = 26,
                   xlab = expression(x[1]), ylab = expression(x[2]),
		   box01 = FALSE, \dots)
}
\arguments{
  \item{x}{a \code{"\linkS4class{Copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{FUN}{the \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{the number of grid points used in each dimension.  This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{a small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{xlab, ylab}{the x-axis and y-axis labels.}
  \item{xlim, ylim}{the \code{\link{range}} of the x and
    y variables, respectively.}
  \item{box01}{a logical specifying if a faint rectangle should be drawn
    on the boundary of \eqn{[0,1]^2} (often useful for copulas, but
    typically \emph{not} for general multivariate distributions
    (\code{"mvdc"})).}
  \item{\dots}{further arguments for (the default method of)
    \code{\link{contour}()}, e.g., \code{nlevels}, \code{levels}, etc.}
}
\section{Methods}{
  Contour lines are drawn for \code{"\linkS4class{Copula}"} or
  \code{"\linkS4class{mvdc}"} objects, see \code{x} in the
  \emph{Arguments} section.
}
%\author{Marius Hofert, Ivan Kojadinovic and Jun Yan}
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
contour(r, levels = seq(1, max(r$z), by=2), lwd=1.5)
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
