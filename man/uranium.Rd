\name{uranium}
\Rdversion{1.1}
\alias{uranium}
\docType{data}
\title{Uranium Exploration Dataset of Cook & Johnson (1986)}
\description{
  These data consist of log concentrations of 7 chemical elements in 655
  water samples collected near Grand Junction, CO (from the Montrose
  quad-rangle of Western Colorado).  Concentrations were
  measured for the following elements: Uranium (U), Lithium (Li), Cobalt
  (Co), Potassium (K), Cesium (Cs), Scandium (Sc), And Titanium (Ti).
}
\usage{data(uranium)}
\format{
  A data frame with 655 observations of the following 7 variables:
  \describe{
    \item{\code{U}}{(numeric) log concentration of Uranium.}
    \item{\code{Li}}{(numeric) log concentration of Lithium.}
    \item{\code{Co}}{(numeric) log concentration of Colbalt.}
    \item{\code{K}}{(numeric) log concentration of Potassium.}
    \item{\code{Cs}}{(numeric) log concentration of Cesium.}
    \item{\code{Sc}}{(numeric) log concentration of Scandum.}
    \item{\code{Ti}}{(numeric) log concentration of Titanium.}
  }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Cook, R. D.  and  Johnson, M. E. (1986)
  Generalized BurrParetologistic distributions with applications to a
  uranium exploration data set.
  \emph{Technometrics} \bold{28}, 123--131.
}
\examples{
data(uranium)
}
\keyword{datasets}
