\name{interval}
\alias{interval}
\title{Construct simple "interval" object}
\description{
  Easy construction of an object of class \code{\linkS4class{interval}},
  using typical mathematical notation.
}
\usage{
interval(ch)
}
\arguments{
  \item{ch}{a character string specifying the interval.}
}
\value{
  an \code{\linkS4class{interval}} object.
}
\author{Martin Maechler}
\seealso{
  the \code{\linkS4class{interval}} class documentation,
  notably its reference to more sophisticated interval classes available
  for \R.
}
\examples{
interval("[0, 1)")

## Two ways to specify open interval borders:
identical(interval("]-1,1["),
          interval("(-1,1)"))

## infinite :
interval("[0, Inf)")

## arithmetic with scalars works:
4 + 2 * interval("[0, 1.5)") # ->  [4, 7)

## str() to look at internals:
str( interval("[1.2, 7]") )
}
\keyword{arith}
\keyword{utilities}
