\name{setTheta}
\title{Specify the Parameter(s) of a Copula}
\alias{setTheta}
\alias{setTheta,acopula,ANY-method}
\alias{setTheta,copula,ANY-method}
\alias{setTheta,xcopula,ANY-method}
\alias{setTheta,outer_nacopula,numeric-method}
\alias{setTheta,ellipCopula,ANY-method}
\alias{setTheta,khoudrajiCopula,ANY-method}
\alias{setTheta,mixCopula,ANY-method}
\description{
  Set or change the parameter \eqn{\theta}{theta} (\code{theta}) of
  a copula.  The name \sQuote{theta} has been from its use in (nested)
  Archimedean copulas, where \code{x} is of class
  \code{"\linkS4class{acopula}"} or
  \code{"\linkS4class{outer_nacopula}"}.  This is used for constructing
  copula models with specified parameter, as, for example, in
  \code{\link{onacopula}()}, or also \code{\link{gofCopula}}.
}
\usage{
setTheta(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)

\S4method{setTheta}{acopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)
\S4method{setTheta}{copula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)
\S4method{setTheta}{xcopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)
\S4method{setTheta}{outer_nacopula,numeric}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)
\S4method{setTheta}{khoudrajiCopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE, \dots)

\S4method{setTheta}{mixCopula,ANY}(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE,
         treat.negative = c("set.0", "warn.set0", "stop"), \dots)
}
\arguments{
  \item{x}{an \R object of class \code{\linkS4class{Copula}}, i.e., any
    copula from package \pkg{copula}.}
  \item{value}{parameter value or vector, \code{\link{numeric}} or
    \code{\link{NA}} (when \code{na.ok} is true), corresponding to the
    \emph{\dQuote{free}} parameters.}
  \item{\dots}{further arguments for methods.}
  \item{na.ok}{logical indicating if \code{\link{NA}} values are ok for
    \code{theta}.}
  \item{noCheck}{logical indicating if parameter constraint checks should
    be skipped.}
  \item{freeOnly}{logical indicating that only non-fixed aka
    \dQuote{free} parameters are to be set.  If true as by default,
    \code{setTheta()} modifies only the free parameters of the copula;
    see also \code{\link{fixParam}}.}
  \item{treat.negative}{a \code{\link{character}} string indicating how
    negative mixture weights should be handled.  If not \code{"stop"}
    which produces an error via \code{\link{stop}}, negative mixture
    weights are replaced by zero.}
}
\value{
  an \R object of the same class as \code{x}, with the main parameter (vector)
  (often called \code{theta}) set to \code{value}.
%% FIXME: be more explicit about outer_nacopula  case
}
\seealso{
  the \dQuote{inverse} function, a \dQuote{getter} method, is
  \code{\link{getTheta}()}.
}
%\author{Martin Maechler}
\examples{
myC <- setTheta(copClayton, 0.5)
myC
## Frank copula with Kendall's tau = 0.8 :
(myF.8 <- setTheta(copFrank, iTau(copFrank, tau = 0.8)))

# negative theta is ok for dim = 2 :
myF <- setTheta(copFrank, -2.5, noCheck=TRUE)
myF@tau(myF@theta) # -0.262

myT <- setTheta(tCopula(df.fixed=TRUE), 0.7)
stopifnot(all.equal(myT, tCopula(0.7, df.fixed=TRUE),
                    check.environment=FALSE, tolerance=0))

(myT2 <- setTheta(tCopula(dim=3, df.fixed=TRUE), 0.7))
## Setting 'rho' and 'df'  --- for default df.fixed=FALSE :
(myT3 <- setTheta(tCopula(dim=3), c(0.7, 4)))
}
\keyword{manip}
