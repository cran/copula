\name{gasoil}
\title{Daily Crude Oil and Natural Gas Prices from 2003 to 2006}
\alias{gasoil}
\docType{data}
\description{
  Three years of daily prices (from July 2003 to July 2006) of crude oil
  and natural gas. These data should be very close to those analysed in
  Gr\enc{é}{e}goire, Genest and Gendron (2008).
}
\usage{data(gasoil, package="copula")}
\format{
  A data frame of 762 daily prices from 2003 to 2006.
  \describe{
    \item{\code{date}}{date (of class \code{\link{Date}}).}
    \item{\code{oil}}{daily price of crude oil}
    \item{\code{gas}}{daily price of natural gas}
  }
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Gr\enc{é}{e}goire, V., Genest, C., and Gendron, M. (2008)
  Using copulas to model price dependence in energy markets.
  \emph{Energy Risk} \bold{5}(5), 58--64.
}
\examples{
data(gasoil)
## Log Scaled  Oil & Gas Prices :
lattice :: xyplot(oil + gas ~ date, data = gasoil, auto.key=TRUE,
                  type = c("l","r"),
                  scales=list(y = list(log = TRUE), equispaced.log = FALSE))
}
\keyword{datasets}
