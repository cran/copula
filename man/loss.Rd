\name{loss}
\Rdversion{1.1}
\alias{loss}
\docType{data}
\title{LOSS and ALAE insurance data}
\description{
Indemnity payment and allocated loss adjustment expense from an insurance
company.
}
\usage{data(loss)}
\format{
  A data frame with 1500 observations of the following 4 variables:
  \describe{
    \item{\code{loss}}{a numeric vector of loss amount up to the \code{limit}.}
    \item{\code{alae}}{a numeric vector of the corresponding allocated loss
  adjustment expense.}
    \item{\code{limit}}{a numeric vector of limit (-99 means no limit).}
    \item{\code{censored}}{1 means censored (limit reached) and 0 otherwise.}
  }
}
\references{
Frees, E. and Valdez, E. (1998). Understanding relationships using
copulas. \emph{North American Actuarial Journal} \bold{2}, 1--25.
}
\examples{
data(loss)
}
\keyword{datasets}
