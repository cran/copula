\name{archmCopula}
\alias{archmCopula}
\alias{claytonCopula}
\alias{frankCopula}
\alias{gumbelCopula}
\alias{amhCopula}
\alias{joeCopula}
\title{Construction of Archimedean Copula Class Object}
\description{
  Constructs an Archimedean copula class object with its corresponding
  parameter and dimension.
}
\usage{
archmCopula(family, param = NA_real_, dim = 2, ...)

claytonCopula(param = NA_real_, dim = 2,
          use.indepC = c("message", "TRUE", "FALSE"))
frankCopula(param = NA_real_, dim = 2,
          use.indepC = c("message", "TRUE", "FALSE"))
gumbelCopula(param = NA_real_, dim = 2,
          use.indepC = c("message", "TRUE", "FALSE"))
amhCopula(param = NA_real_, dim = 2,
          use.indepC = c("message", "TRUE", "FALSE"))
joeCopula(param = NA_real_, dim = 2,
          use.indepC = c("message", "TRUE", "FALSE"))
}
\arguments{
  \item{family}{a character string specifying the family of an
    Archimedean copula.  Currently supported families are "clayton",
    "frank", "amh", "gumbel", and "joe".
  }
  \item{param}{number (\code{\link{numeric}}) specifying the copula parameter.}
  \item{dim}{the dimension of the copula.}
  \item{\dots}{further arguments, passed to the individual creator
    functions (\code{claytonCopula()}, etc).}
  \item{use.indepC}{a string specifying if the independence copula
    \code{\link{indepCopula}}, should be returned in the case where the
    parameter \eqn{\theta}, \code{param}, is at the boundary or limit
    case where the corresponding Archimedean copula is the independence
    copula.
    The default does return \code{indepCopula()} with a message, using
    \code{"TRUE"} does it without a message.  This makes the resulting
    object more useful typically, but does not return a formal
    Archimedean copula of the desired family, something needed e.g., for
    fitting purposes, where you'd use \code{use.indepC="FALSE"}.}
}
\value{
  An Archimedean copula object of class \code{"claytonCopula"},
  \code{"frankCopula"}, \code{"gumbelCopula"}, \code{"amhCopula"}, or \code{"joeCopula"}.
}
\details{
  \code{archmCopula()} is a wrapper for \code{claytonCopula()},
  \code{frankCopula()}, \code{gumbelCopula()}, \code{amhCopula()} and
  \code{joeCopula}.

%% FIXME
  For the Ali-Mikhail-Haq copula family (\code{"amhCopula"}), only the
  bivariate case is available.

  The maximum dimension for which the expression of the pdf is available
  is 6 for the Clayton, Gumbel and Frank families. The cdf expression is
  always available.

  % The maximum dimension for which \code{"dCopula"} can be evaluated is
  % 10 for the Clayton and Gumbel families, and 6 for the Frank family.
  % They are also the maximum dimensions for which maximum likelihood
  % estimation can be done.

  It is now recommended to work with the
  \code{\linkS4class{acopula}}-classed Archimedean copulas, as there is
  no restriction on the dimension there.
}
\references{
  R.B. Nelsen (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{ellipCopula}}, \code{\link{evCopula}}.

  %% -> ./copFamilies.Rd
  \code{\linkS4class{acopula}}-classed Archimedean copulas, such as
  \code{\link{copClayton}}, \code{\link{copGumbel}}, etc, notably for
  mathematical definitions including the meaning of \code{param}.
}
\examples{
clayton.cop <- claytonCopula(2, dim = 3)
## scatterplot3d(rCopula(1000, clayton.cop))

frank.3D  <- frankCopula(dim = 3)# with NA parameters
frank.cop <- frankCopula(3)
persp(frank.cop, dCopula)

gumbel.cop <- archmCopula("gumbel", 5)
contour(gumbel.cop, dCopula)

amh.cop <- amhCopula(0.5)
u. <- as.matrix(expand.grid(u=(0:10)/10, v=(0:10)/10, KEEP.OUT.ATTRS=FALSE))
du <- dCopula(u., amh.cop)
stopifnot(is.finite(du) | apply(u. == 0, 1,any)| apply(u. == 1, 1,any))

## A 7-dim Frank copula
frank.cop <- frankCopula(3, dim = 7)
x <- rCopula(5, frank.cop)
## dCopula now *does* work:
dCopula(x, frank.cop)

## A 7-dim Gumbel copula
gumbel.cop <- gumbelCopula(2, dim = 7)
dCopula(x, gumbel.cop)

## A 12-dim Joe copula
joe.cop <- joeCopula(iTau(joeCopula(), 0.5), dim = 12)
dCopula(x, joe.cop)
}
\keyword{distribution}
\keyword{multivariate}
