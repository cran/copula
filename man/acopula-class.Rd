\name{acopula-class}
\Rdversion{1.1}
\title{Class "acopula" of Archimedean Copula Families}
\docType{class}
\alias{acopula-class}
\alias{acopula}
\alias{initialize,acopula-method}
\alias{show,acopula-method}
\description{
  This class \code{"acopula"} of Archimedean Copula Families
  is mainly used for providing objects of known Archimedean families
  with all related functions.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("acopula", ...)}.
  For several well-known \ifelse{latex}{Archi-medean}{Archimedean}
  copula families, the package \pkg{copula} already provides such
  family objects.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{A string (class \code{"character"}) describing
      the copula family, for example, "AMH" (or simply "A"), "Clayton" ("C"),
      "Frank" ("F"), "Gumbel" ("G"), or "Joe" ("J").}
    \item{\code{theta}:}{Parameter value, a \code{\link{numeric}}, where
      \code{NA} means \dQuote{unspecified}.}
    \item{\code{psi}, \code{iPsi}:}{The (Archimedean) generator
      \eqn{\psi}{psi} (with \eqn{\psi}{psi}(t)=\eqn{\exp}{exp}(-t) being the
      generator of the independence copula) and its inverse
      (\code{\link{function}}). \code{iPsi} has an optional argument
      \code{log} which, if \code{TRUE} returns the logarithm of inverse of the
      generator.}
    \item{\code{absdPsi}:}{A \code{\link{function}} which computes the absolute
      value of the derivative of the generator \eqn{\psi}{psi} for the given
      parameter \code{theta} and of the given \code{degree} (defaults to 1).
      Note that there is no informational loss by computing the absolute value
      since the derivatives alternate in sign (the generator derivative
      is simply (-1)^\code{degree}\eqn{*}\code{absdPsi}).  The
      number \code{n.MC} denotes the sample size for a Monte Carlo evaluation
      approach.  If \code{n.MC} is zero (the default), the
      generator derivatives are evaluated with their exact formulas.  The
      optional parameter \code{log} (defaults to \code{FALSE}) indicates whether
      or not the logarithmic value is returned.}
    \item{\code{absdiPsi}:}{a \code{\link{function}} computing the
      absolute value of the derivative of the generator inverse
      (\code{iPsi()}) for the given parameter \code{theta}.  The optional
      parameter \code{log} (defaults to \code{FALSE}) indicates whether
      the logarithm of the absolute value of the first derivative of
      \code{iPsi()} is returned.}
    \item{\code{dDiag}:}{a \code{\link{function}} computing the density of the
      diagonal of the Archimedean copula at \code{u} with parameter
      \code{theta}.  The parameter \code{log} is as described before.}
    \item{\code{dacopula}:}{a \code{\link{function}} computing the density
      of the Archimedean copula at \code{u} with parameter \code{theta}.  The
      meanings of the parameters \code{n.MC} and \code{log} are as
      described before.}
    \item{\code{score}:}{a \code{\link{function}} computing the
      \emph{derivative} of the density with respect to the parameter
      \eqn{\theta}{theta}.}
    \item{\code{uscore}:}{a \code{\link{function}} computing the
      \emph{derivative} of the density with respect to the each of the arguments.}
    \item{\code{paraInterval}:}{Either \code{\link{NULL}} or an object of class
      \code{"\linkS4class{interval}"}, which is typically obtained from a call
      such as \code{\link{interval}("[a,b)")}.}
    \item{\code{paraConstr}:}{A function of \code{theta} returning
      \code{TRUE} if and only if \code{theta} is a valid parameter value.  Note
      that \code{paraConstr} is built automatically from the interval,
      whenever the \code{paraInterval} slot is valid.
      \code{"\linkS4class{interval}"}.}
    \item{\code{nestConstr}:}{A \code{\link{function}}, which returns
      \code{TRUE} if and only if the two provided parameters \code{theta0} and
      \code{theta1} satisfy the sufficient nesting condition for this family.}
    \item{\code{V0}:}{A \code{\link{function}} which samples \code{n} random
      variates from the distribution \eqn{F} with Laplace-Stieltjes transform
      \eqn{\psi}{psi} and parameter \code{theta}.}
    \item{\code{dV0}:}{A \code{\link{function}} which computes either the
      probability mass function or the probability density function (depending
      on the Archimedean family) of the distribution function whose
      Laplace-Stieltjes transform equals the generator \eqn{\psi}{psi} at the
      argument \code{x} (possibly a vector) for the given parameter
      \code{theta}.  An optional argument \code{log} indicates whether
      the logarithm of the mass or density is computed (defaults to
    	\code{FALSE}).}
    \item{\code{V01}:}{A \code{\link{function}} which gets a vector of
      realizations of \code{V0}, two parameters \code{theta0} and \code{theta1}
      which satisfy the sufficient nesting condition, and which
      returns a vector of the same length as \code{V0} with random variates from
      the distribution function \eqn{F_{01}}{F01} with Laplace-Stieltjes
      transform \eqn{\psi_{01}}{psi01} (see \code{dV01}) and parameters
      \eqn{\theta_0=\,}{theta0=}\code{theta0},
      \eqn{\theta_1=\,}{theta1=}\code{theta1}.}
    \item{\code{dV01}:}{Similar to \code{dV0} with the difference being that
      this \code{\link{function}} computes the probability mass or density
      function for the Laplace-Stieltjes transform
      \deqn{\psi_{01}(t;V_0) = \exp(-V_0\psi_0^{-1}(\psi_1(t))),
      }{psi01(t;V0) = exp(-V0*psi0^{-1}(psi1(t))),}
      corresponding to the distribution function \eqn{F_{01}}{F01}.

      Arguments are the evaluation point(s) \code{x}, the value(s) \code{V0},
      and the parameters \code{theta0} and \code{theta1}.  As for
      \code{dV0}, the optional argument \code{log} can be specified
      (defaults to \code{FALSE}). Note that if \code{x} is a vector,
      \code{V0} must either have length one (in which case \code{V0} is
      the same for every component of \code{x}) or \code{V0} must be of
      the same length as \code{x} (in which case the components of
      \code{V0} correspond to the ones of \code{x}).}
    \item{\code{tau}, \code{iTau}:}{Compute Kendall's tau of the bivariate
      Archimedean copula with generator \eqn{\psi}{psi} as a
      \code{\link{function}} of \code{theta}, respectively, \code{theta} as a
      function of Kendall's tau.}
    \item{\code{lambdaL}, \code{lambdaU}, \code{lambdaLInv},
      \code{lambdaUInv}:}{Compute the lower (upper) tail-dependence coefficient
      of the bivariate Archimedean copula with generator \eqn{\psi}{psi} as a
      \code{\link{function}} of \code{theta}, respectively, \code{theta} as a
      function of the lower (upper) tail-dependence coefficient.}
  }
  For more details about Archimedean families, corresponding distributions and
  properties, see the references.
}
\section{Methods}{
  \describe{
    \item{initialize}{\code{signature(.Object = "acopula")}: is used to
      automatically construct the function slot \code{paraConstr}, when
      the \code{paraInterval} is provided (typically via
      \code{\link{interval}()}).}

    \item{show}{\code{signature("acopula")}: compact overview of the copula.}
  }
}
%\author{Martin Maechler, Marius Hofert}
\seealso{
  Specific provided copula family objects, for example,
  \code{\link{copAMH}}, \code{\link{copClayton}}, \code{\link{copFrank}},
  \code{\link{copGumbel}}, \code{\link{copJoe}}.
  \cr To access these, you may also use \code{\link{getAcop}}.

  A \emph{nested} Archimedean copula \emph{without}
  child copulas (see class \code{"\linkS4class{nacopula}"})
  is a proper \ifelse{latex}{Archime-dean}{Archimedean} copula, and hence,
  \code{\link{onacopula}()} can be used to construct a specific
  parametrized Archimedean copula; see the example below.
}
\references{% -> ./copFamilies.Rd
  See those of the families, for example, \code{\link{copGumbel}}.
}
\examples{
## acopula class information
showClass("acopula")

## Information and structure of Clayton copulas
copClayton
str(copClayton)

## What are admissible parameters for Clayton copulas?
copClayton@paraInterval

## Can two Clayton copulas with parameters theta0 and theta1 be nested?
## Case 1: theta0 = 3, theta1 = 2
copClayton@nestConstr(theta0 = 3, theta1 = 2)
## -> FALSE as the sufficient nesting criterion is not fulfilled
## Case 2: theta0 = 2, theta1 = 3
copClayton@nestConstr(theta0 = 2, theta1 = 3) # TRUE

## For more examples, see  help("acopula-families")
}
\keyword{multivariate}
\keyword{classes}
