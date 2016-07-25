\name{pairs2}
\alias{pairs2}
\title{Scatter-Plot Matrix ('pairs') for Copula Distributions with Nice Defaults}
\description{
  A version of \pkg{graphics}' package \code{\link{pairs}()},
  particularly useful for visualizing dependence in multivariate
  (copula) data.
}
\usage{
pairs2(x, labels = NULL, labels.null.lab = "U", \dots)
}
\arguments{
  \item{x}{a numeric \code{\link{matrix}} or an \R object for which
    \code{\link{as.matrix}(x)} returns such a matrix.}
  \item{labels}{the variable names, typically unspecified.}
  \item{labels.null.lab}{the \code{\link{character}} string
    determining the \dQuote{base name} of the variable labels in case
    \code{labels} is \code{NULL} and \code{x} does not have all column names given.}
  \item{\dots}{further arguments, passed to \code{\link{pairs}()}.}
}
\value{\code{\link{invisible}()}
}
%\author{Marius Hofert}
\seealso{
  \code{\link{splom2}()} for a similar function based on
  \code{\link[lattice]{splom}()}.
}
\examples{
## Create a 100 x 7 matrix of random variates from a t distribution
## with four degrees of freedom and plot the generated data
U <- matrix(rt(700, df = 4), ncol = 7)
pairs2(U, pch = ".")
}
\keyword{hplot}
