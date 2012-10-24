\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,normalCopula-method}
\alias{show,tCopula-method}
\alias{show,copula-method}
\alias{show,fitCopula-method}
\alias{show,fitMvdc-method}
\title{Methods for `show' in Package `copula'}
\description{
  Methods for function \code{\link{show}} in package \pkg{copula}.
}
\section{Methods}{
  \describe{
    \item{object = "copula"}{see \code{\link{Copula}}. }

    \item{object = "fitMvdc"}{(see \code{\link{fitMvdc}}): and}
    \item{object = "fitCopula"}{(see \code{\link{fitCopula}}): these call the
      (hidden) \code{\link{print}} method, with its default argument.
      Using \code{\link{print}()} instead, allows to set \code{digits}, e.g.}
    \item{object = "fitMvdc"}{see \code{\link{fitCopula}}. }
}}
%\examples{

%}
\keyword{methods}
\keyword{print}
