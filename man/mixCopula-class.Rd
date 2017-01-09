\name{mixCopula-class}
\title{Class \code{"mixCopula"} of Copula Mixtures}
\alias{mixCopula-class}
\alias{dim,mixCopula-method}
\alias{lambda,mixCopula-method}
\alias{rho,mixCopula-method}
\Rdversion{1.1}
\docType{class}
\description{
  The class \code{"mixCopula"} is the class of all finite mixtures of
  copulas.

  These are given by (a list of) \eqn{m} \dQuote{arbitrary} copulas, and
  their respective \eqn{m} non-negative probability weights.
}
\section{Objects from the Class}{
  Objects are typically created by \code{\link{mixCopula}()}.
}
\section{Slots}{
  \describe{
    \item{\code{w}:}{Object of class \code{"mixWeights"}, basically a
      non-negative \code{\link{numeric}} vector of length, say \eqn{m},
      which sums to one.}
    \item{\code{cops}:}{Object of class \code{"parClist"}, a
      \code{\link{list}} of (parametrized) copulas, \code{"\linkS4class{parCopula}"}.}
  }
}
\note{
  As the probability weights must some to one (\code{1}), which is part
  of the validity (see \code{\link{validObject}}) of an object of class
  \code{"mixWeights"}, the number of \dQuote{free} parameters inherently
  is (at most) one \emph{less} than the number of mixture components
  \eqn{m}.

  Because of that, it does not make sense to fix (see
  \code{\link{fixParam}} or \code{\link{fixedParam<-}})
  all but one of the weights: Either all are fixed, or at least two must
  be free.  Further note, that the definition of free or fixed
  parameters, and the meaning of the methods (for \code{mixCopula}) of
  \code{\link{getTheta}}, \code{\link{setTheta}} and
  \code{\link{fixedParam<-}} will probably change in a next release of
  package \pkg{copula}, where it is planned to use a reparametrization
  better suited for \code{\link{fitCopula}}.
}
\section{Extends}{
  Class \code{"\linkS4class{parCopula}"}, directly.
  Class \code{"\linkS4class{Copula}"}, by class "parCopula", distance 2.
}
\section{Methods}{
  \describe{
    \item{dim}{\code{signature(x = "mixCopula")}: dimension of copula.}
    \item{rho}{\code{signature(x = "mixCopula")}: Spearman's rho of
      copula \code{x}.}
    \item{lambda}{\code{signature(x = "mixCopula")}: lower and upper
      tail dependecies \code{\link{lambda}}, \eqn{(\lambda[L],\lambda[U])},
      of the mixture copula.}
  }
}
\seealso{
  \code{\link{mixCopula}} for creation and examples.
}
\examples{
showClass("mixCopula")
}
\keyword{classes}
