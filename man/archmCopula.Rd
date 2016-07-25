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

  For the mathematical definitions of the respective Archimedean families,
  see \code{\link{copClayton}}.

  For \eqn{d = 2}, i.e. \code{dim = 2}, the AMH, Clayton and Frank copulas
  allow to model negative Kendall's tau (\code{\link{tau}}) behavior via
  negative \eqn{\theta}{theta}, for AMH and Clayton
  \eqn{-1 \le \theta}{-1 <= theta}, and for Frank
  \eqn{-\infty < \theta}{-Inf < theta}.
  The respective boundary cases are
  \describe{
    \item{AMH:}{\eqn{\tau(\theta = -1) = -0.1817258},}
    \item{Clayton:}{\eqn{\tau(\theta = -1) = -1},}
    \item{Frank:}{\eqn{\tau(\theta = -Inf) = -1} (as limit).}
  }

  For the Ali-Mikhail-Haq copula family (\code{"amhCopula"}), only the
  bivariate case is available; however \code{\link{copAMH}} has no such
  limitation.

  The maximum dimension for which the expression of the pdf is available
  is 6 for the Clayton, Gumbel and Frank families. The cdf expression is
  always available.

  % The maximum dimension for which \code{"dCopula"} can be evaluated is
  % 10 for the Clayton and Gumbel families, and 6 for the Frank family.
  % They are also the maximum dimensions for which maximum likelihood
  % estimation can be done.
  For \eqn{d > 2}, i.e. \code{dim > 2}, it is now recommended to work with the
  \code{\linkS4class{acopula}}-classed Archimedean copulas, as there is
  no restriction on the dimension there.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  R.B. Nelsen (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  %% -> ./copFamilies.Rd
  \code{\linkS4class{acopula}}-classed Archimedean copulas, such as
  \code{\link{copClayton}}, \code{\link{copGumbel}}, etc, notably for
  mathematical definitions including the meaning of \code{param}.

  \code{\link{fitCopula}} for fitting such copulas to data.

  \code{\link{ellipCopula}}, \code{\link{evCopula}}.
}
\examples{
clayton.cop <- claytonCopula(2, dim = 3)
## scatterplot3d(rCopula(1000, clayton.cop))

## negative param (= theta) is allowed for dim = 2 :
tau(claytonCopula(-0.5)) ## = -1/3
tauClayton <- Vectorize(function(theta) tau(claytonCopula(theta, dim=2)))
plot(tauClayton, -1, 10, xlab=quote(theta), ylim = c(-1,1), n=1025)
abline(h=-1:1,v=0, col="#11111150", lty=2); axis(1, at=-1)

tauFrank <- Vectorize(function(theta) tau(frankCopula(theta, dim=2)))
plot(tauFrank, -40, 50, xlab=quote(theta), ylim = c(-1,1), n=1025)
abline(h=-1:1,v=0, col="#11111150", lty=2)

## tauAMH() is function in our package
iTau(amhCopula(), -1) # -1 with a range warning
iTau(amhCopula(), (5 - 8*log(2)) / 3) # -1 with a range warning

ic <- frankCopula(0) # independence copula (with a "message")
stopifnot(identical(ic,
   frankCopula(0, use.indepC = "TRUE")))# indep.copula  withOUT message
(fC <- frankCopula(0, use.indepC = "FALSE"))
## a Frank copula which corresponds to the indep.copula (but is not)

frankCopula(dim = 3)# with NA parameters
frank.cop <- frankCopula(3)# dim=2
persp(frank.cop, dCopula)

gumbel.cop <- archmCopula("gumbel", 5)
stopifnot(identical(gumbel.cop, gumbelCopula(5)))
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
gumbelC.7 <- gumbelCopula(2, dim = 7)
dCopula(x, gumbelC.7)

## A 12-dim Joe copula
joe.cop <- joeCopula(iTau(joeCopula(), 0.5), dim = 12)
dCopula(x, joe.cop)
}
\keyword{distribution}
\keyword{multivariate}
