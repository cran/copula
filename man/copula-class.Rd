\name{copula-class}
\title{Mother Classes "Copula", etc of all Copulas in the Package}
\docType{class}
\alias{copula-class}
\alias{Copula-class}
\alias{dimCopula-class}
\alias{parCopula-class}
\alias{xcopula-class}
% methods
\alias{dim,copula-method}
\alias{dim,xcopula-method}
\alias{dim,dimCopula-method}
\description{
  A copula is a multivariate distribution with uniform margins.
  The virtual class \code{"Copula"} is the mother (or \dQuote{super
    class}) of all copula classes in the package \pkg{copula} which
  encompasses classes of the former packages \pkg{nacopula} and \pkg{copula}.

  The virtual class \code{"parCopula"} extends \code{"Copula"} and is
  the super class of all copulas that can be fully \emph{par}ametrized and
  hence fitted to data.  For these, at least the \code{\link{dim}()} method
  must be well defined.

  The virtual class \code{"dimCopula"} extends \code{"Copula"} and has
  an explicit slot \code{dimension}, with corresponding trivial
  \code{\link{dim}()} method.

  The virtual class \code{"copula"} extends bot \code{"dimCopula"} and
  \code{"parCopula"} and is the mother of all copula classes from former
  package \pkg{copula}.  It has set of slots for (the dimension and)
  parameter vector, see below.

  The virtual class \code{"xcopula"} extends \code{"parCopula"} and
  contains a slot \code{copula}; an (\dQuote{actual}) class example are
  the rotated copulas, \code{\linkS4class{rotCopula}}.
}
\section{Objects from the Class}{
  Objects are typically created by are by \code{\link{tCopula}()},
  \code{\link{evCopula}()}, etc.
}
\section{Slots}{
  Class \code{"dimCopula"} and its subclasses, notably
  \code{"copula"}, have a slot
  \describe{ \item{\code{dimension}:}{an \code{"integer"} (of length 1), the
      copula dimension \eqn{d}.} }

  Class \code{"copula"} (and its subclasses) have \emph{additional} slots
  \describe{
    \item{\code{parameters}:}{\code{\link{numeric}} vector of
      parameter values, can be NA (i.e., \code{\link{NA_real_}}).}
    \item{\code{param.names}:}{\code{"character"} vector of parameter
      names (and hence of the same length as \code{parameters}).}
    \item{\code{param.lowbnd}:}{lower bounds for the parameters, of
      class \code{"numeric"}.}
    \item{\code{param.upbnd}:}{upper bounds for the parameters, of class
      \code{"numeric"}.}
    \item{\code{fullname}:}{\bold{deprecated}; object of class \code{"character"}, family names
      of the copula.}
  }
}
\note{
  The \code{"copula"} class is extended by the
  \code{\linkS4class{evCopula}}, \code{\linkS4class{archmCopula}},
  and \code{\linkS4class{ellipCopula}} classes.  Instances of
  such copulas can be created via functions \code{\link{evCopula}},
  \code{\link{archmCopula}} and \code{\link{ellipCopula}}.

  \code{"plackettCopula"} and \code{\linkS4class{fgmCopula}} are
  special types of copulas which do not belong to either one of the
  three classes above.
}
\section{Warning}{
  This implementation is still at the experimental
  stage and is subject to change during the development.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\seealso{
  Help for the (sub)classes
  \code{\linkS4class{archmCopula}},
  \code{\linkS4class{ellipCopula}},
  \code{\linkS4class{evCopula}}, and
  \code{\linkS4class{fgmCopula}}.

  The Archimedean and nested Archimedean classes (from former package
  \pkg{nacopula}), with a more extensive list of slots (partly instead
  of methods),
  \code{\linkS4class{acopula}}, and \code{\linkS4class{nacopula}}.
}
\examples{
hc <- evCopula("husler", 1.25)
dim(hc)
smoothScatter(u <- rCopula(2^11, hc))
lambda   (hc)
tau (hc)
rho(hc)
str(hc)
}
\keyword{classes}
