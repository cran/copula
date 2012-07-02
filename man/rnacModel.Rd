\name{rnacModel}
\alias{rnacModel}
\title{Random nacopula model}
\description{
  Randomly construct a nested Archimedean copula model,
  %% FIXME
}
\usage{
rnacModel(family, d, pr.comp, rtau0 = function() rbeta(1, 2,4),
          order=c("random", "each", "seq"), digits.theta = 2)
}
\arguments{
  \item{family}{the Archimedean family}
  \item{d}{integer >=2; the dimension}
  \item{pr.comp}{probability of a direct component on each level}
  \item{rtau0}{a \code{\link{function}} to generate a (random) tau,
    corresponding to theta0, the outermost theta.}
  \item{order}{string indicating how the component IDs are selected.}
  \item{digits.theta}{integer specifying the number of digits to round
    the theta values.}
}
% \details{
% }
\value{
  an object of \code{\linkS4class{outer_nacopula}}.
}
\author{Martin Maechler, 10 Feb 2012}
\seealso{
  \code{\link{rnacopula}} for generating \eqn{d}-dimensional
  observations from an (outer) \code{\linkS4class{nacopula}}, e.g., from
  the \emph{result} of \code{rnacModel()}.
}
\examples{
## Implicitly tests the function {with validity of outer_nacopula ..}
set.seed(11)
for(i in 1:40) {
  m1 <- rnacModel("Gumbel", d=sample(20:25, 1), pr.comp = 0.3,
		  rtau0 = function() 0.25)
  m2 <- rnacModel("Joe", d=3, pr.comp = 0.1, order="each")
  mC <- rnacModel("Clayton", d=20, pr.comp = 0.3,
		  rtau0 = function() runif(1, 0.1, 0.5))
  mF <- rnacModel("Frank", d=sample(20:25, 1), pr.comp = 0.3, order="seq")
}
}% Ex.
\keyword{distribution}
\keyword{multivariate}

