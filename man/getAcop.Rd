\name{getAcop}
\title{Get "acopula" Family Object by Name}
\alias{getAcop}
\alias{getAname}
\alias{.ac.shortNames}
\alias{.ac.longNames}
\alias{.ac.objNames}
\alias{.ac.classNames}
\description{
  Get one of our "acopula" family objects (see \link{acopula-families}
  by name.

  Named strings for \dQuote{translation} between different names and
  forms of Archimedean copulas.
}
\usage{
getAcop (family, check = TRUE)
getAname(family, objName = FALSE)

.ac.shortNames
.ac.longNames
.ac.objNames
.ac.classNames
}
\arguments{
  \item{family}{either a \code{\link{character}} string, the short or longer
    form of the Archimedean family name (for example, "Clayton" or simply "C";
    see the \code{\link{acopula-families}} documentation), or an
    \code{\linkS4class{acopula}} family object, or an object inheriting
    from class \code{\linkS4class{archmCopula}}.}
  \item{check}{logical indicating whether the class of the return value
    should be checked to be \code{"\linkS4class{acopula}"}.}
  \item{objName}{logical indicating that the \emph{name} of the \R
    object should be returned, instead of the family name, e.g.,
    \code{"copClayton"} instead of \code{"Clayton"}.}
}
\value{\code{getAcop()} returns an \code{"\linkS4class{acopula}"} family
  object, typically one of one of our predefined ones.

  \code{getAname()} returns a \code{\link{character}} string, the name
  of an \code{"\linkS4class{acopula}"} family object.

  \code{.as.longnames} etc are named string constants, useful in
  programming for all our (five) standard Archimedean families.
}
%\author{Martin Maechler}
\seealso{
  Our predefined \code{\link{acopula-families}};
  the class definition \code{"\linkS4class{acopula}"}.
}
\examples{
getAcop("Gumbel")

## different ways of getting the same "acopula" family object:
stopifnot(## Joe (three ways):
          identical(getAcop("J"), getAcop("Joe")),
          identical(getAcop("J"), copJoe),
          ## Frank (yet another two different ways):
          identical(getAcop(frankCopula()), copFrank),
          identical(getAcop("frankCopula"), copFrank))

stopifnot(
  identical(getAname(claytonCopula()), getAname("C")),
  identical(getAname(copClayton), "Clayton"), identical(getAname("J"), "Joe"),
  identical(getAname(amhCopula(), TRUE), "copAMH"),
  identical(getAname(joeCopula(), TRUE), "copJoe")
)

.ac.shortNames
.ac.longNames
.ac.objNames
.ac.classNames
}
\keyword{manip}
