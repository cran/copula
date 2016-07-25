\name{evTestA}
\alias{evTestA}
\title{Bivariate Test of Extreme-Value Dependence Based on Pickands'
  Dependence Function}
\description{
  Test of bivariate extreme-value dependence based on the process
  comparing the empirical copula with a natural nonparametric
  estimator of the unknown copula derived under extreme-value
  dependence.  The test statistics are defined in the third reference.
  Approximate p-values for the test statistics are obtained
  by means of a \emph{multiplier} technique.
}
\usage{
evTestA(x, N = 1000, derivatives = c("An","Cn"))
}
\arguments{
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{number of multiplier iterations to be used to simulate
    realizations of the test statistic under the null hypothesis.}
  \item{derivatives}{string specifying how the derivatives of the unknown
    copula are estimated, either \code{"An"} or \code{"Cn"}.  The former
    gives better results for samples smaller than 400 but is slower.}
}
\details{
  More details are available in the third reference.
  See also Genest and Segers (2009) and Remillard and Scaillet (2009).
}
\value{
  Returns a list whose attributes are:
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Genest, C. and Segers, J. (2009). Rank-based inference for bivariate
  extreme-value copulas. \emph{Annals of Statistics}, 37, pages 2990-3022.

  R\enc{é}{e}millard, B. and Scaillet, O. (2009). Testing for equality
  between two copulas. \emph{Journal of Multivariate Analysis}, 100(3),
  pages 377-386.

  Kojadinovic, I. and Yan, J. (2010).
  Nonparametric rank-based tests of bivariate extreme-value dependence.
  \emph{Journal of Multivariate Analysis} \bold{101}, 2234--2249.

  Kojadinovic, I. and Yan, J. (2010).
  Modeling Multivariate Distributions with Continuous Margins Using the
  copula R Package.
  \emph{Journal of Statistical Software} \bold{34}(9), 1--20.
  \url{http://www.jstatsoft.org/v34/i09/}.
}
\note{
  This test was derived under the assumption of continuous margins,
  which implies that ties occur with probability zero.  The
  presence of ties in the data might substantially affect the
  approximate p-value.  One way of dealing with ties was suggested in
  the last reference.
}
\seealso{ \code{\link{evTestK}}, \code{\link{evTestC}},
  \code{\link{evCopula}},
  \code{\link{gofEVCopula}}, \code{\link{An}}.
}
\examples{
## Do these data come from an extreme-value copula?
set.seed(63)
uG <- rCopula(100, gumbelCopula(3))
uC <- rCopula(100, claytonCopula(3))
\donttest{## takes time: 48 seconds on MM's lynne (2012-06)
evTestA(uG)
evTestA(uG, derivatives = "Cn")

evTestA(uC)
}
\dontshow{% test only: very small subset to be fast
set.seed(34)
tG   <- evTestA(uG[1:10,])
tC.A <- evTestA(uC[1:12,])
tC.C <- evTestA(uC[1:20,], derivatives = "Cn")
stopifnot(
 all.equal(tG  $p.value, 27/2002,  tolerance= 1e-15),
 all.equal(tC.C$p.value,   1/182 , tolerance= 1e-15),
 all.equal(tC.C$p.value, tC.A$p.value, tolerance = 1e-15) # in < 10\% of cases
)
}%dont
}
\keyword{htest}
\keyword{multivariate}
