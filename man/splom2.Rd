\name{splom2}
\alias{splom2}
\title{Scatterplot Matrix (splom) with Nice Variable Names}
\description{
  A version of \pkg{lattice}'s \code{\link{splom}} function,
  particularly useful for visualizing multivariate data sampled from
  copulas, notably nested Archimedean ones.

  \bold{Experimental} We may replace the interface entirely,
  for example, to accept an \code{"\linkS4class{outer_nacopula}"}.
}
\usage{
splom2(data, varnames=NULL, Vname="U", xlab="",
       col.mat = NULL, bg.col.mat = NULL, \dots)
}
\arguments{
  \item{data}{numeric matrix or as.matrix(.)able.}
  \item{varnames}{variable names, typically unspecified.}
  \item{Vname}{character string to become the "base name" of the variables.}
  \item{xlab}{x-axis label.}
  \item{col.mat}{matrix of colors for the plot symbols (the default is the
    setting as obtained from
    \code{\link[lattice]{trellis.par.get}("plot.symbol")$col}).}
  \item{bg.col.mat}{matrix of colors for the background (the default is the
    setting as obtained from
    \code{\link[lattice]{trellis.par.get}("background")$col}).}
  \item{\dots}{further arguments, passed to \code{\link[lattice]{splom}()}.}
}
\value{from \code{\link{splom}()}, an \R object of class \code{"trellis"}.
}
\author{Martin Maechler, with a hint from Deepayan Sarkar; based on
  ideas from Marius Hofert.}
\examples{
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
U7 <- matrix(rt(700, 4), 100, 7)
G <- splom2(U7)
G
}
\keyword{hplot}
