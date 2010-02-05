\name{uranium}
\Rdversion{1.1}
\alias{uranium}
\docType{data}
\title{
Uranium exploration dataset of Cook & Johnson (1986)
}
\description{
These data consist of 655 chemical analyses from water samples collected
from the Montrose quad-rangle of Western Colorado. Concentrations were
measured for the following elements: uranium (U), Lithium (Li), cobalt
(Co), potassium (K), cesium (Cs), scandium (Sc), and titanium (Ti).
}
\usage{data(uranium)}
\format{
  A data frame with 655 observations of the following 7 variables:
  \describe{
    \item{\code{U}}{a numeric vector of uranium.}
    \item{\code{Li}}{a numeric vector of Lithium.}
    \item{\code{Co}}{a numeric vector of colbalt.}
    \item{\code{K}}{a numeric vector of potassium.}
    \item{\code{Cs}}{a numeric vector of cesium.}
    \item{\code{Sc}}{a numeric vector of scandum.}
    \item{\code{Ti}}{a numeric vector of titanium.}
  }
}
\references{
R. D. Cook & M. E. Johnson (1986). Generalized BurrParetologistic
distributions with applications to a uranium exploration data
set. \emph{Technometrics}, 28:123--131.
}
\examples{
data(uranium)
}
\keyword{datasets}
