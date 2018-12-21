\name{fhCopula-class}
\title{Class "fhCopula" of Fr\enc{é}{e}chet-Hoeffding Bound Copulas}
\docType{class}
\alias{fhCopula-class}
\alias{lowfhCopula-class}
\alias{upfhCopula-class}
\description{
  Fr\enc{é}{e}chet-Hoeffding bound copula class.
}
\section{Objects from the Class}{
  Created by calls of the form \code{new("fhCopula", ...)} or rather
  typically by \code{\link{fhCopula}()}, \code{lowfhCopula()}, or
  \code{upfhCopula()}.  Actual (sub) classes are the lower and upper
  Fr\enc{é}{e}chet-Hoeffding bound copulas \code{lowfhCopula()}
  (\eqn{W}), and \code{upfhCopula()} (\eqn{M}).
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{inherited from
      super class \code{"\linkS4class{dimCopula}"}.}
    \item{\code{exprdist}:}{an \code{\link{expression}} of length two,
      named \code{"cdf"} with the \R expression of the CDF, and
      \code{"pdf"} which is empty as the PDF does not exist (everywhere).}
  }
}
%\author{Marius Hofert}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{evCopula}}.

  The class constructors are \code{\link{fhCopula}()},
  \code{lowfhCopula()}, and \code{upfhCopula()}.  See there for
  examples.
}
\keyword{classes}
