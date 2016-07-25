\name{allComp}
\alias{allComp}
\title{All Components of a (Inner or Outer) Nested Archimedean Copula}
\description{
  Given the nested Archimedean copula \code{x}, return an integer vector
  of the \emph{indices} of all components of the corresponding
  \code{\linkS4class{outer_nacopula}} which are components of \code{x},
  either direct components or components of possible child copulas.  This
  is typically only used by programmers investigating the exact nesting
  structure.

  For an \code{\linkS4class{outer_nacopula}} object
  \code{x}, \code{allComp(x)} must be the same as
  \code{1:\link{dim}(x)}, whereas its \dQuote{inner} component copulas
  will each contain a \emph{subset} of those indices only.
}
\usage{
allComp(x)
}
\arguments{
  \item{x}{an \R object inheriting from class \code{\linkS4class{nacopula}}.}
}
\value{
  An \code{\link{integer}} vector of indices \eqn{j} of all components
  \eqn{u_j} as described in the description above.
}
%\author{Martin Maechler}
\examples{
 C3 <- onacopula("AMH", C(0.7135, 1, C(0.943, 2:3)))
 allComp(C3) # components are 1:3
 allComp(C3@childCops[[1]]) # for the child, only  (2, 3)
}
\keyword{manip}
\keyword{utilities}
