\name{setTheta}
\alias{setTheta}
\alias{setTheta,acopula,ANY-method}
\alias{setTheta,outer_nacopula,numeric-method}
\alias{setTheta,copula,ANY-method}
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
\S4method{setTheta}{outer_nacopula,numeric}(x, value, na.ok = TRUE, noCheck = FALSE, \dots)
}
\arguments{
  \item{x}{acopula}
  \item{value}{parameter value or vector, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true.)}
  \item{\dots}{further arguments for methods.}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
  \item{noCheck}{logical indicating if parameter contraint checks should
    be skipped.}
}
\value{
  an \R object of the same class as \code{x}, with the main parameter (vector)
  (often called \code{theta}) set to \code{value}.
%% FIXME: be more explicit about outer_nacopula  case
}
\details{%% FIXME
  Buglet: Not yet working sensibly for \code{tCopulat(df.fixed=FALSE)}.
}
\author{Martin Maechler}
\examples{
myC <- setTheta(copClayton, 0.5)
myC

# negative theta is ok for dim = 2 :
myF <- setTheta(copFrank, -2.5, noCheck=TRUE)
myF@tau(myF@theta) # -0.262

myT <- setTheta(tCopula(), 0.7)
stopifnot(all.equal(myT, tCopula(0.7), tol=0))
}
\keyword{manip}
