\name{show-methods}
\title{Methods for 'show()' in Package 'copula'}
\docType{methods}
\alias{show-methods}
\alias{show,normalCopula-method}
\alias{show,tCopula-method}
\alias{show,parCopula-method}
\alias{show,fitCopula-method}
\alias{show,fitMvdc-method}
\description{
  Methods for function \code{\link{show}} in package \pkg{copula}.
}
\section{Methods}{
  \describe{
    \item{object = "parCopula"}{see \code{\link{Copula}}. }

    \item{object = "fitMvdc"}{(see \code{\link{fitMvdc}}): and}
    \item{object = "fitCopula"}{(see \code{\link{fitCopula}}): these call the
      (hidden) \code{\link{print}} method, with its default argument.
      Using \code{\link{print}()} instead, allows to set \code{digits}, e.g.}
    \item{object = "fitMvdc"}{see \code{\link{fitCopula}}. }
}}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
%\examples{
%}
\keyword{methods}
\keyword{print}
