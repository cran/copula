\name{khoudrajiCopula-class}
\title{Class \code{"khoudrajiCopula"} and its Subclasses}
\docType{class}
\alias{asymCopula-class}
\alias{asym2Copula-class}
\alias{khoudrajiCopula-class}
\alias{khoudrajiBivCopula-class}
\alias{khoudrajiExplicitCopula-class}
%
\alias{dim,khoudrajiCopula-method}

\description{
  The \emph{virtual} class \code{"asymCopula"} of (conceptually) all
  asymmetric copulas and its \sQuote{subclass} \code{"asym2copula"} of
  those which are constructed from two copulas.

  More specifically, the class \code{"khoudrajiCopula"} and its two
  subclasses \code{"khoudrajiBivCopula"} and \code{"khoudrajiExplicitCopula"}
  represent copulas constructed using Khoudraji's device from two copulas of
  the same dimension; see \code{\link{khoudrajiCopula}()} for more details.
}
\section{Objects from the Class}{
  Objects are typically created via \code{\link{khoudrajiCopula}(...)}.
}
\section{Slots}{
  As these classes extend \code{"\linkS4class{copula}"}, they have all
  its slots: \code{dimension}, \code{parameters}, \code{param.names},
  \code{param.lowbnd}, \code{param.upbnd}, and \code{fullname}.
  The classes \code{"khoudrajiCopula"} and \code{"khoudrajiBivCopula"}
  have the extra slots
  \describe{
    \item{\code{copula1}:}{object of class \code{"\linkS4class{copula}"}.}
    \item{\code{copula2}:}{second object of class  \code{"\linkS4class{copula}"}.}
  }
  In addition to these, the class \code{"khoudrajiExplicitCopula"} has the slots
  \describe{
    \item{\code{exprdist}:}{an \code{\link{expression}}, ...}% FIXME
    \item{\code{derExprs1}:}{an \code{\link{expression}} of length \eqn{d}, ...}% FIXME
    \item{\code{derExprs2}:}{an \code{\link{expression}} of length \eqn{d}, ...}% FIXME
  }
}
\section{Methods}{
  When possible, methods are defined at the \code{"khoudrajiCopula"}
  class level.  The implementation of method \code{\link{dCopula}} for
  instance is however not possible at that level.  In addition, it
  differs for \code{"khoudrajiBivCopula"} and
  \code{"khoudrajiExplicitCopula"} classes.
}
\references{
  Genest, C., Ghoudi, K., and Rivest, L.-P. (1998), Discussion of
  "Understanding relationships using copulas", by Frees, E., and
  Valdez, E., \emph{North American Actuarial Journal} \bold{3},
  143--149.

  Khoudraji, A. (1995), Contributions à l'\enc{é}{e}tude des copules et
  \enc{à}{a}la mod\enc{é}{e}lisation des valeurs extr\enc{ê}{e}mes
  bivari\enc{é}{e}es, \emph{PhD thesis, Université Laval}, Québec,
  Canada.

  Liebscher, E. (2008), Construction of asymmetric multivariate copulas,
  \emph{Journal of Multivariate Analysis} \bold{99}, 2234--2250.
}
%% \author{Ivan Kojadinovic}
\seealso{ \code{\link{khoudrajiCopula}()}}
\examples{
showClass("khoudrajiCopula")# two subclasses

## all methods applicable to these subclasses:
(meths <- sapply(names(getClass("khoudrajiCopula")@subclasses),
                 function(CL) methods(class = CL),
                 simplify=FALSE))
}
\keyword{classes}
