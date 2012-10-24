\name{nesdepth}
\alias{nesdepth}
\title{Nesting Depth of a Nested Archimedean Copula ("nacopula")}
\description{
  Compute the nesting depth of a nested Archimedean copula which is the
  length of the longest branch in the tree representation of the copula,
  and hence at least one.
}
\usage{
nesdepth(x)
}
\arguments{
  \item{x}{object of class \code{"\linkS4class{nacopula}"}.}
}
\value{an integer, the nesting depth of the nested Archimedean copula.
  An (unnested) Archimedean copula has depth \code{1}.
}
\seealso{
  \code{\link[=nacopula-class]{dim}} of nacopulas.
}
\examples{
F2 <- onacopula("F", C(1.9, 1, C(4.5, c(2,3))))
F2
F3 <- onacopula("Clayton", C(1.5, 3:1,
                             C(2.5, 4:5,
                               C(15, 9:6))))
nesdepth(F2) # 2
nesdepth(F3) # 3
\dontshow{
stopifnot(identical(nesdepth(F2), 2L),
          identical(nesdepth(F3), 3L),
          identical(nesdepth(onacopula("Gumbel", C(1.5, 3:1))), 1L))
}
}
\keyword{utilities}
