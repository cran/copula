\name{printNacopula}
\alias{printNacopula}
\alias{show,nacopula-method}
\title{Print Compact Overview of a Nested Archimedean Copula ("nacopula")}
\description{
  Print a compact overview of a nested Archimedean copula, that is, an
  object of class \code{"\linkS4class{nacopula}"}.
  Calling \code{printNacopula} explicitly allows to customize the
  printing behavior.  Otherwise, the \code{\link{show}()} method calls
  \code{printNacopula} with default arguments only.
}
\usage{
printNacopula(x, labelKids=NA, deltaInd=, indent.str="",
	              digits=getOption("digits"),
	              width=getOption("width"), ...)
}
\arguments{
  \item{x}{an \R object of class \code{\linkS4class{nacopula}}.}
  \item{labelKids}{logical specifying if child copulas should be labeled;
    If \code{NA} (as per default), on each level, children are labeled
    only if they are not only-child.}
  \item{deltaInd}{by how much should each child be indented \emph{more}
    than its parent? (non-negative integer).  The default is three
    with \code{labelKids} being the default or \code{TRUE}, otherwise it is
    five (for \code{labelKids=FALSE}).}
  \item{indent.str}{a \code{\link{character}} string specifying the
    indentation, that is, the string that should be \emph{prepended} on the
    first line of output, and determine the amount of blanks for the
    remaining lines.}
  \item{digits, width}{number of significant digits, and desired print
    width; see \code{\link{print.default}}.}
  \item{\dots}{potentially further arguments, passed to methods.}
}
\value{
  invisibly, \code{x}.
}
\author{Martin Maechler}
\examples{
C8 <- onacopula("F", C(1.9, 1,
                       list(K1 = C(5.7, c(2,5)),
                            abc= C(5.0, c(3,4,6),
                                   list(L2 = C(11.5, 7:8))))))
C8 # -> printNacopula(C8)
printNacopula(C8, delta=10)
printNacopula(C8, labelKids=TRUE)
}
\keyword{utilities}
