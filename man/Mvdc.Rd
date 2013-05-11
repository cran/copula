\name{Mvdc}
\alias{Mvdc}
\alias{mvdc}
\alias{dMvdc}
\alias{pMvdc}
\alias{rMvdc}
\alias{dmvdc}%deprecated
\alias{pmvdc}%deprecated
\alias{rmvdc}%deprecated
\title{Multivariate Distributions Constructed from Copulas}
\description{
  Density, distribution function, and random generator for a
  multivariate distribution via copula.
}
\usage{
mvdc(copula, margins, paramMargins, marginsIdentical = FALSE,
     check = TRUE, fixupNames = TRUE)
dMvdc(x, mvdc, log=FALSE)
pMvdc(x, mvdc)
rMvdc(n, mvdc)
}
\arguments{
  \item{copula}{ an object of \code{"\linkS4class{copula}"}. }
  \item{margins}{ a character vector specifying all the marginal
    distributions.  See details below.}
  \item{paramMargins}{a \code{\link{list}} whose each component is a
    list (or numeric vectors) of named components, giving the parameter
    values of the marginal distributions.  See details below.}
  \item{marginsIdentical}{ logical variable restricting the marginal
    distributions to be identical. }
  \item{check}{logical indicating to apply quick checks about existence of
    \code{margins} \dQuote{p*} and \dQuote{d*} functions.}
  \item{fixupNames}{logical indicating if the parameters of the margins
    should get automatic names (from \code{\link{formals}(p<mar_i>)}).}
  \item{mvdc}{a \code{"\linkS4class{mvdc}"} object. }
  \item{x}{ a vector of the copula dimension or a matrix with number of
    columns being the copula dimension, giving the coordinates of the
    points where the density or distribution function needs to be
    evaluated.}
  \item{log}{logical indicating if the \code{\link{log}} density should
    be returned.}
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
  marginal distributions. For example, the list
  \deqn{\code{paramMargins = list(list(mean = 0, sd = 2), list(rate = 2))}}
  can be used to specify that the first margin is normal with mean 0 and
  standard deviation 2, and the second margin is exponential with rate 2.
}
\value{
  \code{mvdc()} constructs an object of class \code{"\linkS4class{mvdc}"}.
  \code{dMvdc()} gives the density, \code{pMvdc()} gives the cumulative
  distribution function, and \code{rMvdc()} generates random variates.
}
% \references{
%   Yan (2007). Enjoy the joy of copulas. \emph{Journal of Statistical
%     Software}, 21(4): 1--21.
% }
\seealso{
  \code{\link{ellipCopula}},
  \code{\link{archmCopula}};
  the classes \code{\linkS4class{mvdc}} and \code{\linkS4class{copula}}.
}
\examples{
## construct a bivariate distribution whose marginals
## are normal and exponential respectively, coupled
## together via a normal copula
mv.NE <- mvdc(normalCopula(0.75), c("norm", "exp"),
              list(list(mean = 0, sd =2), list(rate = 2)))
dim(mv.NE)
mv.NE  # using its print() / show() method

# Generate (bivariate) random numbers from that, and visualize
x.samp <- rMvdc(100, mv.NE)
plot(x.samp)
dMvdc(x.samp, mv.NE)
pMvdc(x.samp, mv.NE)
persp  (mv.NE, dMvdc, xlim = c(-4, 4), ylim=c(0, 1), main = "dMvdc(mv.NE)")
persp  (mv.NE, pMvdc, xlim = c(-4, 4), ylim=c(0, 1), main = "pMvdc(mv.NE)")
contour(mv.NE, dMvdc, xlim = c(-4, 4), ylim=c(0, 1))
}
\keyword{distribution}
\keyword{multivariate}
