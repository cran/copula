\name{persp-methods}
\docType{methods}
\alias{persp-methods}
\alias{persp,Copula-method}
\alias{persp,mvdc-method}
\title{Methods for Function `persp' in Package `copula'}
\description{
  Methods for function \code{\link{persp}} to draw perspective plots (of
  two dimensional distributions from package \pkg{copula}).
}
\usage{
\S4method{persp}{Copula}(x, FUN, n.grid = 26, delta = 0,
      xlab = "u1", ylab = "u2",
      zlab = deparse(substitute(FUN))[1], zlim = NULL,
      theta = -30, phi = 30, expand = 0.618,
      ticktype = "detail", \dots)
\S4method{persp}{mvdc}(x, FUN, xlim, ylim, n.grid = 26,
      xlab = "x1", ylab = "x2", zlab = deparse(substitute(FUN))[1],
      theta = -30, phi = 30, expand = 0.618,
      ticktype = "detail", \dots)
}
\arguments{
  \item{x}{a \code{"\linkS4class{Copula}"} or a \code{"\linkS4class{mvdc}"} object.}
  \item{FUN}{the \code{\link{function}} to be plotted; typically
    \code{\link{dCopula}} or \code{\link{pCopula}}.}
  \item{n.grid}{the number of grid points used in each dimension.  This
    can be a vector of length two, giving the number of grid points used
    in x- and y-direction, respectively; the function \code{FUN} will be
    evaluated on the corresponding (x,y)-grid.}
  \item{delta}{A small number in \eqn{[0, \frac 1 2)}{[0, 1/2)}
    influencing the evaluation boundaries. The x- and y- vectors will
    have the range \code{[0+delta, 1-delta]}, the default being \code{[0,1]}.}
  \item{xlim, ylim}{The \code{\link{range}} of the x and
    y variables, respectively.}
  \item{xlab, ylab, zlab, zlim, theta, phi, expand, ticktype, \dots}{Arguments
    for (the default method of) \code{\link{persp}()}, the ones
    enumerated here all with \emph{different} defaults than there.}
}
\section{Methods}{
  Perspective plots for both \code{"\linkS4class{copula}"} or
  \code{"\linkS4class{mvdc}"} objects, see \code{x} in the
  \emph{Arguments} section.
}
\value{
  \code{\link{invisible}}; a list with the following components:
  \item{x, y}{The numeric vectors, as passed to \code{\link{persp.default}}.}
  \item{z}{The matrix of evaluated \code{FUN} values on the grid as
    passed to \code{\link{persp.default}}.}
  \item{persp}{the \eqn{4 \times 4}{4 x 4} transformation matrix
    returned by \code{\link{persp.default}}.}
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{
  The \code{\link{contour-methods}} for drawing contour lines of the
  same functions.
}
\examples{
persp(claytonCopula(2),   pCopula, main = "CDF of claytonCopula(2)")
persp(  frankCopula(1.5), dCopula, main = "Density of frankCopula(1.5)")
persp(  frankCopula(1.5), dCopula, main = "c_[frank(1.5)](.)", zlim = c(0,2))

## Examples with negative tau:
(th1 <- iTau(amhCopula(), -0.1))
persp(amhCopula(th1), dCopula)
persp(amhCopula(th1), pCopula, ticktype = "simple") # no axis ticks
persp(  frankCopula(iTau(  frankCopula(), -0.1)), dCopula)
persp(claytonCopula(iTau(claytonCopula(), -0.1)), dCopula)
##
cCop.2 <- function(u, copula, ...) cCopula(u, copula, ...)[,2]
persp(    amhCopula(iTau(    amhCopula(), -0.1)), cCop.2, main="cCop(AMH...)[,2]")
persp(  frankCopula(iTau(  frankCopula(), -0.1)), cCop.2, main="cCop(frankC)[,2]")
## and  Clayton  also looks "the same" ...

## MVDC Examples ------------------------------------
mvNN <- mvdc(gumbelCopula(3), c("norm", "norm"),
          list(list(mean = 0, sd = 1), list(mean = 1)))
persp(mvNN, dMvdc, xlim=c(-2, 2), ylim=c(-1, 3), main = "Density")
persp(mvNN, pMvdc, xlim=c(-2, 2), ylim=c(-1, 3), main = "Cumulative Distr.")
}
\keyword{methods}
\keyword{hplot}
