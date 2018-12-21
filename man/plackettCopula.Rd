\name{plackettCopula}
\alias{plackettCopula}
\title{Construction of a Plackett Copula}
\description{
  Constructs a Plackett copula (class \code{"plackettCopula}"}) with its
  corresponding parameter.
}
\usage{
plackettCopula(param)
}
\arguments{
  \item{param}{a number (numeric vector of length one) specifying the
    \emph{non negative} parameter.}
}
\value{
  A Plackett copula object of class \code{"\linkS4class{plackettCopula}"}.
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Plackett, R. L. (1965). A Class of Bivariate Distributions.
  \emph{Journal of the American Statistical Association} \bold{60}, 516--522.
}
\seealso{
  The \code{"\linkS4class{plackettCopula}"} class;
  \code{\link{ellipCopula}}, \code{\link{archmCopula}}.
}
\examples{
plackett.cop <- plackettCopula(param=2)
lambda(plackett.cop) # 0 0 : no tail dependencies
## For really large param values (here, 1e20 and Inf are equivalent):
set.seed(1); Xe20 <- rCopula(5000, plackettCopula(1e20))
set.seed(1); Xinf <- rCopula(5000, plackettCopula(Inf))
stopifnot(all.equal(Xe20, Xinf))
}
\keyword{distribution}
\keyword{multivariate}
