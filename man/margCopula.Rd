\name{margCopula}
\title{Marginal copula of a Copula With Specified Margins}
\alias{margCopula}
\alias{margCopula,archmCopula,logical-method}
\alias{margCopula,normalCopula,logical-method}
\alias{margCopula,tCopula,logical-method}
\description{
  The marginal copula of a copula \eqn{C(u_1,\dots, u_d)} is simply the
  restriction of \eqn{C} on a subset of the the coordinate (directions) \eqn{u_1,\dots,u_d}.
}
\usage{
margCopula(copula, keep)
}
\arguments{
  \item{copula}{a \code{"\linkS4class{copula}"} (\R object) of
    dimension, \eqn{d}, say.}
  \item{keep}{logical vector (of length \code{d}) indicating which margins to keep.}
}
\details{
  The marginal copula of a copula is needed in practical data analysis when
  one or more of the components of some multivariate observations is
  missing.  For normal/t/Archimedean copulas, the marginal copulas can be
  easily obtained.  For a general copula, this may not be an easy problem.

  The current implementation only supports normal/t/Archimedean
  copulas. \code{margCopula} is generic function with methods for the
  different copula classes.
}
\value{
  The marginal copula of the specified margin(s).
}
%\author{Jun Yan}
\examples{
tc <- tCopula(8:2 / 10, dim = 8, dispstr = "toep")
margCopula(tc, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))

nc <- normalCopula(.8, dim = 8, dispstr = "ar1")
mnc <- margCopula(nc, c(TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
mnc7 <- margCopula(nc, (1:8) != 1)
stopifnot(dim(nc) == 8, dim(mnc) == 4, dim(mnc7) == 7)

gc <- gumbelCopula(2, dim = 8)
margCopula(gc, c(TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE))
}
\keyword{models}
\keyword{multivariate}
