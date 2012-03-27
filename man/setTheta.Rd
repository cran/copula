\name{setTheta}
\alias{setTheta}
\title{Specify the Parameter of an Archimedean Copula}
\description{
  Set or change the parameter \eqn{\theta}{theta} (\code{theta}) of
  an Archimedean copula \code{x} of class \code{"\linkS4class{acopula}"}.  This
  is used for constructing copula models with specified parameter, as, for
  example, in \code{\link{onacopula}()}.
}
\usage{
setTheta(x, value, na.ok = TRUE)
}
\arguments{
  \item{x}{acopula}
  \item{value}{parameter value, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true.)}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
}
\value{
  an \code{"\linkS4class{acopula}"} as \code{x}, just with the slot
  \code{theta} set to \code{value}.
}
\author{Martin Maechler}
\examples{
myC <- setTheta(copClayton, 0.5)
myC
}
\keyword{manip}
