\name{fitCopula-class}
\title{Classes of Fitted Multivariate Models: Copula, Mvdc}
\docType{class}
\alias{fittedMV-class}
\alias{fitCopula-class}
\alias{fitMvdc-class}
\alias{summaryFitCopula-class}
\alias{summaryFitMvdc-class}
% methods --- see also ./show-methods.Rd, and ./fitMvdc.Rd (S3 meth.s)
\alias{summary,fitCopula-method}
\alias{summary,fitMvdc-method}
\description{
  Classes and summary methods related to copula model fitting.
}
\section{Objects from the Class}{
  Objects can be created by calls to \code{\link{fitCopula}} or
  \code{\link{fitMvdc}}, respectively or to their \code{summary} methods.
}
\section{Slots}{
  The \dQuote{mother class}, \code{"fittedMV"} has the slots
  \describe{
    \item{\code{estimate}:}{\code{numeric}, the estimated parameters.}
    \item{\code{var.est}:}{\code{numeric}, variance matrix estimate of
      the parameter estimator. See note below.}
    \item{\code{loglik}:}{\code{numeric}, log likelihood evaluated at
      the maximizer.}
    \item{\code{nsample}:}{\code{numeric}, integer representing the
      sample size.}
    \item{\code{method}:}{\code{character}, method of estimation.}
    \item{\code{fitting.stats}:}{a \code{\link{list}}, currently
      containing the numeric \code{convergence} code from
      \code{\link{optim}}, the \code{counts}, \code{message}, and all
      the \code{control} arguments explicitly passed to
      \code{\link{optim}}.}
  }
  In addition, the \code{"fitCopula"} class has a slot
  \describe{
    \item{\code{copula}:}{the \emph{fitted} copula, of class
      \code{"\linkS4class{copula}"}.}
  }
  whereas the \code{"fitMvdc"} has
  \describe{
    \item{\code{mvdc}:}{the \emph{fitted} distribution, of class
      \code{"\linkS4class{mvdc}"}.}
  }
}
\section{Extends}{
  Classes \code{"fitCopula"} and \code{"fitMvdc"} extend class
  \code{"fittedMV"}, directly.
}
\section{Methods}{
  \describe{
    \item{summary}{\code{signature(object = "fitMvdc")}: ... }
    \item{summary}{\code{signature(object = "fitCopula")}: ... }
  }
  Further, there are S3 methods (class \code{"fittedMV"}) for
  \code{\link{coef}()}, \code{\link{vcov}()} and \code{\link{logLik}()},
  see \code{\link{fitMvdc}}.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Genest, C., Ghoudi, K., and Rivest, L.-P. (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika} \bold{82}, 543--552.
}
\keyword{classes}
