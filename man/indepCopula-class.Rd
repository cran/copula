\name{indepCopula-class}
\docType{class}
\alias{indepCopula-class}
\title{Class "indepCopula" }
\description{Independence copula class.}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("indepCopula",
    ...)} or by function \code{"indepCopula"}. Such objects can be
  useful as special cases of parametric copulas, bypassing
  copula-specific computations such as distribution, density, and sampler.
}
\section{Slots}{
  \describe{
    \item{\code{exprdist}:}{Object of class \code{"expression"}:
      expressions of the cdf and pdf of the copula. These expressions
      are used in function 'pcopula' and 'dcopula'. }
    \item{\code{dimension}:}{Object of class \code{"numeric"}, dimension
    of the copula. }
    \item{\code{parameters}:}{Object of class \code{"numeric"},
      parameter values. }
    \item{\code{param.names}:}{Object of class \code{"character"},
      parameter names. }
    \item{\code{param.lowbnd}:}{Object of class \code{"numeric"},
      parameter lower bounds. }
    \item{\code{param.upbnd}:}{Object of class \code{"numeric"},
      parameter upper bounds. }
    \item{\code{message}:}{Object of class \code{"character"}, family names
      of the copula. }
  }
}
\section{Methods}{
  \describe{
    \item{Afun}{\code{signature(copula = "indepCopula")}: ... }
    \item{dcopula}{\code{signature(copula = "indepCopula")}: ... }
    \item{pcopula}{\code{signature(copula = "indepCopula")}: ... }
    \item{rcopula}{\code{signature(copula = "indepCopula")}: ... }
  }
}
% \section{Extends}{
%   Class \code{"indepCopula"} does not extend class \code{"copula"}.
%   It is simpler than \code{"copula"}.
% }
\seealso{
  \code{\link{indepCopula}}, 
  \code{\link{copula-class}}.
}
\keyword{classes}
