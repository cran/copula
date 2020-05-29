\name{getIniParam}
\title{Get Initial Parameter Estimate for Copula}
\alias{getIniParam}
\alias{getIniParam,mixCopula-method}
\alias{getIniParam,parCopula-method}
\description{
  A (S4) generic function and methods providing a typically cheap method to get
  valid parameters for the \code{copula}, given the \code{data}.  This is
  used, e.g., in \code{\link{fitCopula}()} when \code{start} is not
  specified.
}
\usage{
getIniParam(copula, data, default, named = TRUE, \dots)
}
\arguments{
  \item{copula}{a \code{"\linkS4class{copula}"} object.}
  \item{data}{an \eqn{n\times d}{n x d}-matrix of data to which the copula should be fitted.}
  \item{default}{a parameter vector of correct length, to be used when no
    method is available or the method does \dQuote{not work}.}
  \item{named}{\code{\link{logical}} indicating if the result should have
    \code{\link{names}}.}
  \item{\dots}{optional further arguments to underlying methods.}
}

\section{Methods}{
  \describe{
    \item{\code{signature(copula = "parCopula")}}{Close to a \emph{default}
      method (as class \code{"\linkS4class{parCopula}"} contains most
      copulas), currently mainly trying to use a version of
      \code{\link{fitCopula}(*, method = "itau")} (itself based on moment
      matching \code{\link{iTau}()}).}
    \item{\code{signature(copula = "mixCopula")}}{
      a relatively simple method, for the copula parameters, trying
      \code{getInitParam(cop[[k]])} for each component, and using equal weights \code{w[k]}.}
  }
}
%% \details{
%% }
\value{
  a \code{\link{numeric}} vector of correct length, say \code{param}, which
  should e.g., \dQuote{work} in \code{\link{loglikCopula}(param, u = data, copula)}.
}
%% \author{Martin Maechler}
\seealso{
  \code{\link{getTheta}()} gets such a vector \emph{from} a \code{copula} object;
  \code{\link{fitCopula}}, \code{\link{loglikCopula}}.
}
\examples{
 # TODO !
}
\keyword{models}
