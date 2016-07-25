\name{plackettCopula}
\alias{plackettCopula}
\title{Construction of a Plackett Copula Class Object}
\description{
  Constructs a Plackett copula class object with its corresponding
  parameter.
}
\usage{
plackettCopula(param)
}
\arguments{
  \item{param}{ a numeric vector specifying the parameter values. }
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
  \code{\link{ellipCopula}}, \code{\link{archmCopula}}.
}
\examples{
plackett.cop <- plackettCopula(param=2)
}
\keyword{distribution}
\keyword{multivariate}
