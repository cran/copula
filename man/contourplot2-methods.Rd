\name{contourplot2-methods}
\title{Contour Plot  Methods 'contourplot2' in Package 'copula'}
\alias{contourplot2}
\alias{contourplot2-methods}
\alias{contourplot2,matrix-method}
\alias{contourplot2,data.frame-method}
\alias{contourplot2,Copula-method}
\alias{contourplot2,mvdc-method}
\docType{methods}
\description{
  Methods for \code{contourplot2()}, a version of
  \code{\link[lattice]{contourplot}()} from \pkg{lattice}, to draw
  contour plots of two dimensional distributions from package \pkg{copula}.
}
\usage{
\S4method{contourplot2}{matrix}(x, aspect = 1,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      xlab = NULL, ylab = NULL,
      cuts = 16, labels = !region,
      scales = list(alternating = c(1,1), tck = c(1,0)),
      region = TRUE, \dots,
      col.regions = gray(seq(0.4, 1, length.out = 100)))
\S4method{contourplot2}{data.frame}(x, aspect = 1,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      xlab = NULL, ylab = NULL,
      cuts = 16, labels = !region,
      scales = list(alternating = c(1,1), tck = c(1,0)),
      region = TRUE, \dots,
      col.regions = gray(seq(0.4, 1, length.out = 100)))
\S4method{contourplot2}{Copula}(x, FUN, n.grid = 26, delta = 0,
      xlim = 0:1, ylim = 0:1,
      xlab = quote(u[1]), ylab = quote(u[2]), \dots)
\S4method{contourplot2}{mvdc}(x, FUN, n.grid = 26, xlim, ylim,
      xlab = quote(x[1]), ylab = quote(x[2]), \dots)
}
\arguments{
  \item{x}{a \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{aspect}{the aspect ratio.}
  \item{xlim, ylim}{the x- and y-axis limits.}
  \item{xlab, ylab}{the x- and y-axis labels. If at least one is
    \code{NULL}, useful \code{xlab} and \code{ylab} are determined
    automatically; the behavior depends on the class of \code{x}.}
  \item{cuts}{the number of levels; see \code{\link{contourplot}()}.}
  \item{labels}{a \code{\link{logical}} indicating whether the contour
    lines are labeled; see \code{\link{contourplot}()}.}
  \item{scales}{a \code{\link{list}} determining how the axes are drawn;
    see \code{\link{contourplot}()}.}
  \item{region}{a \code{\link{logical}} indicating whether regions
    between contour lines should be filled as in a level plot; see
    \code{\link{contourplot}()}.}
  \item{col.regions}{the colors of the regions if \code{region = TRUE};
    see \code{\link{contourplot}()}.}
  \item{FUN}{the \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{the number of grid points used in each dimension.  This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{a small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{contourplot}()}.}
}
\section{Methods}{
  Contourplot plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{contourplot}()}.
}
%\author{Marius Hofert}
\seealso{
  The \code{\link{contour-methods}} for drawing perspective plots via base graphics.

  The \pkg{lattice}-based \code{\link{wireframe2-methods}} for functions,
  and \code{\link{cloud2-methods}} and \code{\link{splom2-methods}} for data.
}
\examples{
## For 'matrix' objects
## The Frechet--Hoeffding bounds W and M
n.grid <- 26
u <- seq(0, 1, length.out = n.grid)
grid <- expand.grid("u[1]" = u, "u[2]" = u)
W <- function(u) pmax(0, rowSums(u)-1) # lower bound W
M <- function(u) apply(u, 1, min) # upper bound M
x.W <- cbind(grid, "W(u[1],u[2])" = W(grid)) # evaluate W on 'grid'
x.M <- cbind(grid, "M(u[1],u[2])" = M(grid)) # evaluate M on 'grid'
contourplot2(x.W) # contour plot of W
contourplot2(x.M) # contour plot of M

## For 'Copula' objects
cop <- frankCopula(-4)
contourplot2(cop, pCopula) # the copula
contourplot2(cop, pCopula, xlab = "x", ylab = "y") # adjusting the labels
contourplot2(cop, dCopula) # the density

## For 'mvdc' objects
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
             list(list(mean = 0, sd = 1), list(mean = 1)))
xl <- c(-2, 2)
yl <- c(-1, 3)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, contour = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, region = FALSE, labels = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl, region = FALSE)
contourplot2(mvNN, FUN = dMvdc, xlim = xl, ylim = yl,
             col.regions = colorRampPalette(c("royalblue3", "maroon3"), space="Lab"))
}
\keyword{methods}
\keyword{hplot}
