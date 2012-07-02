\name{nacopula-class}
\Rdversion{1.1}
\docType{class}
\alias{outer_nacopula-class}
\alias{nacopula-class}
\alias{dim,nacopula-method}
\title{Class "nacopula" of nested Archimedean copulas}
\description{
  Class of nested Archimedean Copulas, \code{"nacopula"},
  and its \emph{specification} \code{"outer_nacopula"}
  differ only by the validation method, which is stricter for the
  outer(most) copula (the root copula).
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("nacopula", ...)},
  which is only intended for experts.  Root copulas are typically constructed by
  \code{\link{onacopula}(.)}.
}
\section{Slots}{
  \describe{
    \item{\code{copula}:}{an object of class
      \code{"\linkS4class{acopula}"}, denoting the top-level Archimedean
      copula of the nested Archimedean copula, that is, the root copula.}
    \item{\code{comp}:}{an \code{\link{integer}} vector (possibly of
      length 0) of indices of components in \code{1:d} which are not
      nested Archimedean copulas.  Here, \eqn{d} denotes the
      \emph{dimension} of the random vectors under consideration; see
      the \code{dim()} method below.}
    \item{\code{childCops}:}{a (possibly empty) \code{\link{list}} of
      further nested Archimedean copulas (child copulas), that is, objects
      of class \code{"nacopula"}.  The \code{"nacopula"} objects
      therefore contain \code{"\linkS4class{acopula}"} objects as
      special cases.}
  }
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "nacopula"):} returns the dimension
      \eqn{d} of the random vector \eqn{U} following \code{x}.}
    \item{show}{\code{signature("nacopula"):} calling
      \code{\link{printNacopula}} for a compact overview of the nested
      Archimedean copula under consideration.}
  }
}
\author{Martin Maechler}
\seealso{
  \code{\link{onacopula}} for building (outer) "nacopula" objects.  For the
  class definition of the \code{copula} component, see
  \code{\linkS4class{acopula}}.
}
\examples{
## nacopula and outer_nacopula class information
showClass("nacopula")
showClass("outer_nacopula")

## Construct a three-dimensional nested Frank copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copFrank@tauInv(.2)
theta1 <- copFrank@tauInv(.5)
C3 <- onacopula("F", C(theta0, 1, C(theta1, c(2,3))))

C3 # displaying it, using show(C3); see help(printNacopula)

## What is the dimension of this copula?
dim(C3)

## What are the indices of direct components of the root copula?
C3@comp

## How does the list of child nested Archimedean copulas look like?
C3@childCops # only one child for this copula, components 2, 3

}
\keyword{multivariate}
\keyword{classes}
