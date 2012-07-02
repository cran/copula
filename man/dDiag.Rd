\name{dDiag}
\alias{dDiag}
\title{Density of the diagonal of (nested) Archimedean copulas}
\description{
  Evaluate the density of the diagonal of a \eqn{d}-dimensional (nested)
  Archimedean copula.  Note that the diagonal of a copula is a cumulative
  distribution function.  Currently, only Archimedean copulas are implemented.
}
\usage{
dDiag(u, cop, log=FALSE)
}
\arguments{
  \item{u}{a numeric vector of evaluation points.}
  \item{cop}{a (nested) Archimedean copula object of class
    \code{"\linkS4class{outer_nacopula}"}.  This also determines the dimension
 		via the \code{comp} slot}
  \item{log}{logical indicating if the \code{\link{log}} of the density of the
    diagonal should be returned instead of just the diagonal density.}
}
\value{
  A \code{\link{numeric}} vector containing the values of the density of the
  diagonal of the Archimedean copula at \code{u}.}
\author{Martin Maechler, Marius Hofert}
\references{
	Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
  \emph{Estimators for Archimedean copulas in high dimensions: A comparison},
   submitted.
}
\seealso{
  \code{\linkS4class{acopula}} class, \code{\link{dnacopula}}.
}
\examples{
th. <- c(0.1, 0.2, 0.5, 0.8, 1.4, 2., 5.)
curve(dDiag(x, cop=onacopulaL("Clayton", list(th.[1], 1:3))), 0, 1,
      n=1000, ylab="dDiag(x, *)", main="Diagonal densities of Clayton")
abline(h=0, lty=3)
for(j in 2:length(th.))
  curve(dDiag(x, cop=onacopulaL("Clayton", list(th.[j], 1:3))), add=TRUE,
	     col=j, n=1000)
legend("topleft", do.call(expression, lapply(th., function(th)
                                 substitute(theta == TH, list(TH=th)))),
       lty = 1, col=seq_along(th.), bty="n")
}
\keyword{distribution}
