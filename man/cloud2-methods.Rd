\name{cloud2-methods}
\alias{cloud2}
\alias{cloud2-methods}
\alias{cloud2,matrix-method}
\alias{cloud2,data.frame-method}
\alias{cloud2,Copula-method}
\alias{cloud2,mvdc-method}
\docType{methods}
\title{Cloud Plot Methods ('cloud2') in Package 'copula'}
\description{
  Function and Methods \code{cloud2()} to draw (\CRANpkg{lattice})
  \code{\link{cloud}} plots of two-dimensional distributions from
  package \pkg{copula}.
}
\usage{
\S4method{cloud2}{matrix}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{cloud2}{data.frame}(x,
      xlim = range(x[,1], finite = TRUE),
      ylim = range(x[,2], finite = TRUE),
      zlim = range(x[,3], finite = TRUE),
      xlab = NULL, ylab = NULL, zlab = NULL,
      scales = list(arrows = FALSE, col = "black"),
      par.settings = standard.theme(color = FALSE), \dots)
\S4method{cloud2}{Copula}(x, n,
      xlim = 0:1, ylim = 0:1, zlim = 0:1,
      xlab = quote(U[1]), ylab = quote(U[2]), zlab = quote(U[3]), \dots)
\S4method{cloud2}{mvdc}(x, n,
      xlim = NULL, ylim = NULL, zlim = NULL,
      xlab = quote(X[1]), ylab = quote(X[2]), zlab = quote(X[3]), \dots)
}
\arguments{
  \item{x}{a \code{"\linkS4class{matrix}"},
    \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
    a \code{"\linkS4class{mvdc}"} object.}
  \item{xlim, ylim, zlim}{the x-, y- and z-axis limits.}
  \item{xlab, ylab, zlab}{the x-, y- and z-axis labels.}
  \item{scales}{a \code{\link{list}} determining how the axes are drawn;
    see \code{\link{cloud}()}.}
  \item{par.settings}{see \code{\link{cloud}()}.}
  \item{n}{when \code{x} is not matrix-like: The sample size of the
    random sample drawn from \code{x}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{cloud}()}.}
}
\section{Methods}{
  Cloud plots for objects of class \code{"\linkS4class{matrix}"} ,
  \code{"\linkS4class{data.frame}"}, \code{"\linkS4class{Copula}"} or
  \code{"\linkS4class{mvdc}"}.
}
\value{
  An object of class \dQuote{trellis} as returned by \code{\link[lattice]{cloud}()}.
}
%\author{Marius Hofert; tweaks by MM}
\examples{
## For 'matrix' objects
cop <- gumbelCopula(2, dim = 3)
n <- 1000
set.seed(271)
U <- rCopula(n, copula = cop)
cloud2(U, xlab = quote(U[1]), ylab = quote(U[2]), zlab = quote(U[3]))

## For 'Copula' objects
set.seed(271)
cloud2(cop, n = n) # same as above

## For 'mvdc' objects
mvNN <- mvdc(cop, c("norm", "norm", "exp"),
             list(list(mean = 0, sd = 1), list(mean = 1), list(rate = 2)))
cloud2(mvNN, n = n)
}
\keyword{methods}
\keyword{hplot}
