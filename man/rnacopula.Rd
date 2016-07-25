\name{rnacopula}
\alias{rnacopula}
\title{Sampling Nested Archimedean Copulas}
\description{
  Random number generation for nested Archimedean copulas (of class
  \code{\linkS4class{outer_nacopula}}, specifically), aka
  \emph{sampling} nested Archimedean copulas will generate \code{n}
  random vectors of dimension \eqn{d} (= \code{dim(x)}).
}
\usage{
rnacopula(n, copula, x, ...)
}
\arguments{
  \item{n}{integer specifying the sample size, that is, the number of
    copula-distributed random vectors \eqn{\mathbf{U}_i}{U_i}, to be generated.}
  \item{copula}{an \R object of \code{\link{class}}
    \code{"\linkS4class{outer_nacopula}"}, typically from
    \code{\link{onacopula}()}.}
  \item{x}{only for back compatibility: former name of \code{copula} argument.}
  \item{\dots}{possibly further arguments for the given copula family.}
}
\details{
  The generation happens by calling \code{\link{rnchild}()} on
  each child copula (which itself recursively descends the tree implied
  by the nested Archimedean structure).  The algorithm is based on a
  mixture representation of the generic distribution functions
  \eqn{F_{0}}{F0} and \eqn{F_{01}}{F01} and is presented in
  McNeil (2008) and Hofert (2011a).  Details about how to efficiently sample the
  distribution functions \eqn{F_{0}}{F0} and \eqn{F_{01}}{F01} can be found in
  Hofert (2010), Hofert (2012), and Hofert and \enc{Mächler}{Maechler} (2011).
}
\value{
  \code{\link{numeric}} matrix containing the generated vectors of
  random variates from the nested Archimedean copula object \code{copula}.
}
%\author{Marius Hofert and Martin Maechler}
\references{
  McNeil, A. J. (2008).
  Sampling nested Archimedean copulas.
  \emph{Journal of Statistical Computation and Simulation} \bold{78}, 6,
  567--581.

  Hofert, M. (2010).
  Efficiently sampling nested Archimedean copulas.
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.

  Hofert, M. (2011a).
  A stochastic representation and sampling algorithm for nested
  Archimedean copulas.
  \emph{Journal of Statistical Computation and Simulation}, in press.

  Hofert, M. (2012).
  Sampling exponentially tilted stable distributions.
  \emph{ACM Transactions on Modeling and Computer Simulation} \bold{22},
  1 (3rd article).

  Hofert, M. and \enc{Mächler}{Maechler}, M. (2011).
  Nested Archimedean Copulas Meet R: The nacopula Package.
  \emph{Journal of Statistical Software} \bold{39}, 9, 1--20.
}
\seealso{
  \code{\link{rnchild}}; classes \code{"\linkS4class{nacopula}"} and
  \code{"\linkS4class{outer_nacopula}"}; see also \code{\link{onacopula}()}.
  \code{\link{rnacModel}} creates random nacopula \emph{models}, i.e.,
  the input \code{copula} for \code{rnacopula(n, copula)}.

  Further, those of the Archimedean families, for example,
  \code{\link{copGumbel}}.
}
\examples{
## Construct a three-dimensional nested Clayton copula with parameters
## chosen such that the Kendall's tau of the respective bivariate margins
## are 0.2 and 0.5 :
C3 <- onacopula("C", C(copClayton@iTau(0.2), 1,
                       C(copClayton@iTau(0.5), c(2,3))))
C3
\dontshow{
stopifnot(nrow(rnacopula(1, C3)) == 1,
          nrow(rnacopula(0, C3)) == 0)
}
## Sample n vectors of random variates from this copula.  This involves
## sampling exponentially tilted stable distributions
n <- 1000
U <- rnacopula(n, C3)

## Plot the drawn vectors of random variates
splom2(U)
}
\keyword{distribution}
