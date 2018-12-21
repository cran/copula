\name{moCopula-class}
\docType{class}
\alias{moCopula-class}
\title{Class "moCopula" of Marshall-Olkin Copulas}
\description{
  The Marshall-Olkin copula class.

  The (2-dimensional) "MO" copula with parameter
  \eqn{\bold\theta\in [0,1]^2}{theta in [0,1]^2} is (i.e., its CDF is)
  \deqn{C(u_1, u_2) = min(u_1 * u_2^(1 - \theta_2), u_1^(1 - \theta_1) * u_2).%
  }{C(u1, u2) = min(u1 * u2^(1 - theta2), u1^(1 - theta1) * u2).}
  Consequently, the density is undefined on a curve (in \eqn{[0,1]^2}),
  namely for the points \eqn{\bold{u}=(u_1,u_2)}{u=(u1,u2)} where two
  expressions in the above \eqn{min(f(u), g(u))} are equal,
  \eqn{f(u)=g(u)}.  It is easy to see that that is equivalent to
  %% \eqn{u_1 * u_2^(1 - \theta_2) =  u_1^(1 - \theta_1) * u_2,%
  %% }{u1 * u2^(1  - theta2) =  u1^(1 - theta1) * u,2} which is equivalent to
  \deqn{u_1^{\theta_1} = u_2^{\theta_2}.}{u1^theta_1 = u2^theta2.}
  %% and it is everywhere else, i.e., "almost everywhere"!
}
\section{Objects from the Class}{
  Objects can be created by \code{new("moCopula", ...)} but are
  typically produced by \code{\link{moCopula}(...)}.
}
\section{Slots}{
  \describe{
    \item{\code{dimension}:}{Numeric (scalar), the dimension of the copula.}
    \item{\code{exprdist}:}{a length two \code{\link{expression}} with
      expressions for the CDF and PDF of the copula.}% FIXME: still used?
    \item{\code{parameters}:}{numeric vector of two parameter values in \eqn{[0,1]}.}
    \item{\code{param.names}:}{\code{"\link{character}"} vector of length two.}
    \item{\code{param.lowbnd}:}{numeric vector of two values in \eqn{[0,1]}.}
    \item{\code{param.upbnd}:}{numeric vector of two values in \eqn{[0,1]}.}
    \item{\code{fullname}:}{(deprecated; do not use!)}
  }
}
\section{Methods}{
  Typical copula methods work, see \code{"\linkS4class{moCopula}"}
  and use \code{methods(class = "moCopula")}.
}
\section{Extends}{
  Class \code{"moCopula"} extends class \code{"\linkS4class{copula}"} directly.
}
%\author{Marius Hofert}
\references{
  Nelsen, R. B. (2006), \emph{An introduction to Copulas}, Springer, New York.
}
\seealso{
  \code{\link{moCopula}} for constructing them;
  \code{\link{copula-class}}.
}
\examples{
moCopula()@exprdist[["cdf"]] # a simple definition

methods(class = "moCopula")
contourplot2(moCopula(c(.1, .8)), pCopula, main= "moCopula((0.1, 0.8))")

Xmo <- rCopula(5000, moCopula(c(.2, .5)))
try( # gives an error, as there is no density (!):
loglikCopula(c(.1, .2), Xmo, moCopula())
)

plot(moCopula(c(.9, .2)), n = 10000, xaxs="i", yaxs="i",
     # opaque color (for "density effect"):
     pch = 16, col = adjustcolor("black", 0.3))
}
\keyword{classes}
