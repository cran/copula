\name{empCopula-class}
\title{Class "empCopula" of Empirical Copulas}
\docType{class}
\alias{empCopula-class}
\alias{dim,empCopula-method}% not yet documented
\description{
  Empirical Copula class.
}
\section{Objects from the Class}{
  Created by calls of the form \code{new("empCopula", ...)} or rather
  typically by \code{\link{empCopula}()} based on a matrix
  \code{X} of pseudo-observations.  Smoothing options are available, see
  there.
}
\section{Slots}{
  \describe{
    \item{\code{X}:}{\code{\link{matrix}} of pseudo-observations based
      on which the empirical copula is constructed.}
    \item{\code{smoothing}:}{\code{\link{character}} string determining
      the smoothing method.}
    \item{\code{offset}:}{\code{\link{numeric}} giving the shift in the
      normalizing factor for computing the empirical copula.}% FIXME:
					% length 1 or 2 ?!
    \item{\code{ties.method}:}{a string indicating \code{\link{rank}()}'s
      ties method for computing the empirical copula.}
  }
}
%\author{Marius Hofert}
\seealso{
  The class constructor are \code{\link{empCopula}()}, also for
  examples.
}
\keyword{classes}
