\name{rnchild}
\title{Sampling Child 'nacopula's}
\alias{rnchild}
\description{
  Method for generating vectors of random numbers of nested Archimedean
  copulas which are child copulas.
}
\usage{
rnchild(x, theta0, V0, ...)
}
\arguments{
  \item{x}{an \code{"\linkS4class{nacopula}"} object, typically emerging from an
    \code{"\linkS4class{outer_nacopula}"} object constructed with
    \code{\link{onacopula}()}.}
  \item{theta0}{the parameter (vector) of the parent Archimedean copula
    which contains \code{x} as a child.}
  \item{V0}{a \code{\link{numeric}} vector of realizations of
    \eqn{V_{0}}{V0} following \eqn{F_{0}}{F0} whose length determines the
    number of generated vectors, that is, for each realization
    \eqn{V_{0}}{V0}, a vector of variates from \code{x} is generated.}
  \item{\dots}{possibly further arguments for the given copula family.}
}
\details{
  The generation is done recursively, descending the tree implied by the
  nested Archimedean structure.  The algorithm is based on a mixture
  representation and requires sampling \eqn{V_{01}\sim F_{01}}{V01 ~ F01}
  given random variates \eqn{V_0\sim F_{0}}{V0 ~ F0}.  Calling
  \code{"rnchild"} is only intended for experts.  The typical call of
  this function takes place through \code{\link{rnacopula}()}.
}
\value{
  a list with components
  \item{U}{a \code{\link{numeric}} matrix containing the vector of random
    variates from the child copula.  The number of rows of this matrix
    therefore equals the length of \eqn{V_{0}}{V0} and the number of
    columns corresponds to the dimension of the child copula.}

 \item{indcol}{an \code{\link{integer}} vector of indices of \code{U}
   (the vector following a nested Archimedean copula of which \code{x} is
   a child) whose corresponding components of \code{U} are arguments of
   the nested Archimedean copula \code{x}.}
}
\author{Marius Hofert, Martin Maechler}
\seealso{
  \code{\link{rnacopula}}, also for the references.
  Further, classes \code{"\linkS4class{nacopula}"} and
  \code{"\linkS4class{outer_nacopula}"}; see also \code{\link{onacopula}()}.
}
\examples{
## Construct a three-dimensional nested Clayton copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5.
theta0 <- copClayton@tauInv(.2)
theta1 <- copClayton@tauInv(.5)
C3 <- onacopula("C", C(theta0, 1, C(theta1, c(2,3))))
## Sample n random variates V0 ~ F0 (a Gamma(1/theta0,1) distribution)
n <- 1000
V0 <- copClayton@V0(n, theta0)

## Given these variates V0, sample the child copula, that is, the bivariate
## nested Clayton copula with parameter theta1
U23 <- rnchild(C3@childCops[[1]], theta0, V0)
\dontshow{% also works with 0-length V0:
stopifnot(dim(rnchild(C3@childCops[[1]], theta0, numeric(0))$U) == c(0,2))
}
## Now build the three-dimensional vectors of random variates by hand
U1 <- copClayton@psi(rexp(n)/V0, theta0)
U <- cbind(U1, U23$U)

## Plot the vectors of random variates from the three-dimensional nested
## Clayton copula
splom2(U)
}
\keyword{distribution}
