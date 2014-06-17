\name{persp-methods}
\docType{methods}
\alias{persp-methods}
\alias{persp,copula-method}
\alias{persp,mvdc-method}
\title{Methods for Function `persp' in Package `copula'}
\description{
  Methods for function \code{\link{persp}} to draw perspective plots (of
  two dimensional distributions from package \pkg{copula}.
}
\usage{% want \usage{} since we have many "surprising arguments"
\S4method{persp}{copula}(x, fun,
      n = 51, delta = 0,
      xlab = "x", ylab = "y", zlab = deparse(substitute(fun))[1],
      theta = -30, phi = 30, expand = 0.618,
      ticktype = "detail", \dots)

\S4method{persp}{mvdc}(x, fun,
      xlim, ylim, nx = 51, ny = 51,
      xis = seq(xlim[1], xlim[2], length = nx),
      yis = seq(ylim[1], ylim[2], length = ny),
      xlab = "x", ylab = "y", zlab = deparse(substitute(fun))[1],
      theta = -30, phi = 30, expand = 0.618,
      ticktype = "detail", \dots)
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
    defaulting to zero.  The x- and y- vectors will use range
    \code{[0+delta, 1-delta]}, i.e., \code{[0,1]} by default.}

  %% "mvdc" only :
  \item{xlim, ylim}{(\code{"mvdc"}:) the \code{\link{range}} of the x or
    y variable, respectively.}
  \item{nx,ny}{(\code{"mvdc"}:) the number of points in x- or
    y-direction, respectively.  The function \code{fun} will be
    evaluated on a grid of size \eqn{nx \times ny}{nx x ny}.}
  \item{xis, yis}{(\code{"mvdc"}:) instead of specifying \code{xlim,
      ylim} and \code{nx, ny}, the numeric vectors (of length \code{nx}
    and \code{ny}) may be specified directly.}

  %% persp
  \item{xlab, ylab, zlab, theta, phi, expand, ticktype, \dots}{arguments
    for (the default method of) \code{\link{persp}()}, the ones
    enumerated here all with \emph{different} defaults than there.}
}
\section{Methods}{
  Perspective plots for both \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"} objects, see \code{x} in the
  \emph{Arguments} section.
}
\value{
  invisibly (\code{\link{invisible}}) a list with components
  \item{x,y}{the numeric vectors, as passed to \code{\link{persp.default}}.}
  \item{z}{the \eqn{nx \times ny}{nx x ny} matrix, as passed to
    \code{\link{persp.default}}.}
  \item{persp}{the \eqn{4 \times 4}{4 x 4} transformation matrix
    returned by \code{\link{persp.default}}.}
}
\seealso{
  The \code{\link{contour-methods}} for drawing contour lines of the
  same functions.
}
\examples{
persp(frankCopula(-0.8), dCopula)
persp(claytonCopula(2),  pCopula, main = "CDF of claytonCopula(2)")
## example with negative tau :
(th1 <- iTau(amhCopula(), -0.1))
persp(amhCopula(th1), dCopula)
persp(amhCopula(th1), pCopula, ticktype = "simple")# no axis ticks

mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd =1), list(mean = 1)))
persp(mvNN, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3), main="Density")
persp(mvNN, pMvdc, xlim=c(-2, 2), ylim=c(-1, 3), main="Cumulative Distr.")
}
\keyword{methods}
\keyword{hplot}
