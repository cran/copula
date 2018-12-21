\name{describeCop}
\title{Copula (Short) Description as String}
\alias{describeCop}
\alias{describeCop-methods}
\alias{describeCop,Copula,missing-method}
\alias{describeCop,copula,character-method}
\alias{describeCop,xcopula,ANY-method}
\alias{describeCop,archmCopula,character-method}
\alias{describeCop,ellipCopula,character-method}
\alias{describeCop,fgmCopula,character-method}
\alias{describeCop,moCopula,character-method}
\alias{describeCop,indepCopula,character-method}
\alias{describeCop,khoudrajiCopula,character-method}
\alias{describeCop,mixCopula,character-method}
\alias{describeCop,rotCopula,character-method}
\alias{describeCop,fhCopula,character-method}
\alias{describeCop,empCopula,character-method}
\docType{methods}
\description{
  Describe a \code{\linkS4class{copula}} object, i.e., its basic
  properties as a string.  This is a utility used when
  \code{\link{print}()}ing or \code{\link{plot}()}ting copulas, e.g.,
  after a fitting.
}
\usage{
describeCop(x, kind = c("short", "very short", "long"), prefix = "", ...)
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "archmCopula", kind = "ANY")}}{ .. }
    \item{\code{signature(x = "copula", kind = "character")}}{ .. }
    \item{\code{signature(x = "copula", kind = "missing")}}{ .. }
    \item{\code{signature(x = "ellipCopula", kind = "character")}}{ .. }
    \item{\code{signature(x = "fgmCopula", kind = "ANY")}}{ .. }
    \item{\code{signature(x = "xcopula", kind = "ANY")}}{ .. }
  }
}
\arguments{
  \item{x}{a \code{\linkS4class{copula}} object, or a generalization
    such as  \code{\linkS4class{parCopula}}.}
  \item{kind}{a \code{\link{character}} string specifying the size (or
    \dQuote{complexity} of the copula description desired.}
  \item{prefix}{a string to be prefixed to the returned string, which
    can be useful for indentation in describing extended copulas such
    as Khoudraji copulas.}
  \item{...}{further arguments; unused currently.}
}
%% \details{
%% }
\value{
  a \code{\link{character}} string.
}
\seealso{
  Copula class definition \code{\linkS4class{copula}};
}
\examples{
## FIXME
}
\keyword{utilities}
\keyword{methods}
