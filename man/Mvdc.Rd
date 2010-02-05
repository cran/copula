\name{Mvdc}
\alias{Mvdc}
\alias{mvdc}
\alias{dmvdc}
\alias{pmvdc}
\alias{rmvdc}
\title{Multivariate distributions constructed from copulas}
\description{
  Density, distribution function, and random generator for a
  multivariate distribution via copula.
}
\usage{
mvdc(copula, margins, paramMargins, marginsIdentical = FALSE)
dmvdc(mvdc, x)
pmvdc(mvdc, x)
rmvdc(mvdc, n)
}
\arguments{
  \item{copula}{ an object of \code{"copula"}. }
  \item{margins}{ a character vector specifying all the marginal
    distributions. See details below.}
  \item{paramMargins}{ a list whose each component is a list of named
    components, giving the parameter values of the marginal
    distributions. See details below.}
  \item{marginsIdentical}{ logical variable restricting the marginal
    distributions to be identical. }
  \item{mvdc}{ a \code{"mvdc"} object. }
  \item{x}{ a vector of the copula dimension or a matrix with number of
    columns being the copula dimension, giving the coordinates of the
    points where the density or distribution function needs to be
    evaluated.}
  \item{n}{ number of observations to be generated. }
}
\details{
  The characters in argument \code{margins} are used to construct
  density, distribution, and quantile function names. For
  example, \code{norm} can be used to specify marginal distribution,
  because \code{dnorm}, \code{pnorm}, and \code{qnorm} are all
  available.
  
  A user-defined distribution, for example, \code{fancy}, can be used as
  margin provided that \code{dfancy}, \code{pfancy}, and \code{qfancy} are
  available.
    
  Each component list in argument \code{paramMargins} is a list with
  named components which are used to specify the parameters of the
  marginal distributions. For example,
  \code{paramMargins = list(list(mean = 0, sd = 2), list(rate = 2))}
  can be used to specify that the first margin is normal with mean 0 and
  standard deviation 2, and the second margin is exponential with rate 2.
}
\value{
  \code{"mvdc"} constructs an object of class \code{"mvdc"}.
  \code{"dmvdc"} gives the density, \code{"pmvdc"} gives the distribution
  function, and \code{"rmvdc"} generates random variates.
}
% \references{
%   Yan (2007). Enjoy the joy of copulas. \emph{Journal of Statistical
%     Software}, 21(4): 1--21.
% }
\seealso{
  \code{\link{ellipCopula}},
  \code{\link{archmCopula}},
  \code{\link{mvdc-class}},
  \code{\link{copula-class}}.
}
\examples{
## construct a bivariate distribution whose marginals
## are normal and exponential respectively, coupled
## together via a normal copula
x <- mvdc(normalCopula(0.75), c("norm", "exp"),
          list(list(mean = 0, sd =2), list(rate = 2)))
x.samp <- rmvdc(x, 100)
dmvdc(x, x.samp)
pmvdc(x, x.samp)
persp(x, dmvdc, xlim = c(-4, 4), ylim=c(0, 1))
}
\keyword{ distribution }% at least one, from doc/KEYWORDS
\keyword{ multivariate }% __ONLY ONE__ keyword per line
