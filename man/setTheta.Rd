\name{setTheta}
\alias{setTheta}
\alias{setTheta,acopula,ANY-method}
\alias{setTheta,outer_nacopula,numeric-method}
\alias{setTheta,copula,ANY-method}
\alias{setTheta,ellipCopula,ANY-method}
\title{Specify the Parameter(s) of a Copula}
\description{
%% FIXME:
  Set or change the parameter \eqn{\theta}{theta} (\code{theta}) of
  a copula.  The name \sQuote{theta} has been from its use in (nested)
  Archimedean copulas, where \code{x} is of class
  \code{"\linkS4class{acopula}"} or
  \code{"\linkS4class{outer_nacopula}"}.  This is used for constructing
  copula models with specified parameter, as, for example, in
  \code{\link{onacopula}()}.
}
\usage{
setTheta(x, value, \dots)
\S4method{setTheta}{acopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, \dots)
\S4method{setTheta}{copula,ANY} (x, value, na.ok = TRUE, noCheck = FALSE, \dots)
\S4method{setTheta}{outer_nacopula,numeric}(x, value, na.ok = TRUE, noCheck = FALSE, \dots)
}
\arguments{
  \item{x}{an \R object of class \code{\linkS4class{Copula}}, i.e., any
    copula from package \pkg{copula}.}
  \item{value}{parameter value or vector, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true.)}
  \item{\dots}{further arguments for methods.}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
  \item{noCheck}{logical indicating if parameter constraint checks should
    be skipped.}
}
\value{
  an \R object of the same class as \code{x}, with the main parameter (vector)
  (often called \code{theta}) set to \code{value}.
%% FIXME: be more explicit about outer_nacopula  case
}
%\author{Martin Maechler}
\examples{
myC <- setTheta(copClayton, 0.5)
myC

# negative theta is ok for dim = 2 :
myF <- setTheta(copFrank, -2.5, noCheck=TRUE)
myF@tau(myF@theta) # -0.262

myT <- setTheta(tCopula(df.fixed=TRUE), 0.7)
stopifnot(all.equal(myT, tCopula(0.7, df.fixed=TRUE), tolerance=0))

(myT2 <- setTheta(tCopula(dim=3, df.fixed=TRUE), 0.7))
## Setting 'rho' and 'df'  --- for default df.fixed=FALSE :
(myT3 <- setTheta(tCopula(dim=3), c(0.7, 4)))

}
\keyword{manip}
