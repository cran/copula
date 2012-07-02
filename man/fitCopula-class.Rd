\name{fitCopula-class}
\docType{class}
\alias{fittedMV-class}
\alias{fitCopula-class}
\alias{fitMvdc-class}
\alias{summaryFitCopula-class}
\alias{summaryFitMvdc-class}
\title{Classes of Fitted Multivariate Models: Copula, Mvdc}
\description{
  Classes and summaries related to copula model fitting.
}
\section{Objects from the Class}{
  Objects can be created by calls to \code{\link{fitCopula}},
  \code{\link{fitMvdc}} or to their \code{summary} method.
}
\section{Slots}{
  The \dQuote{mother class}, \code{"fittedMV"} has the slots
  \describe{
    \item{\code{estimate}:}{\code{numeric}, parameter estimate.}
    \item{\code{var.est}:}{\code{numeric}, variance matrix estimate of
      the parameter estimator. See note below.}
    \item{\code{loglik}:}{\code{numeric}, loglikelihood evaluated at
      the maximizer.}
    \item{\code{nsample}:}{\code{numeric}, integer representing the
      sample size.}
    \item{\code{method}:}{\code{character}, method of estimation.}
    \item{\code{convergence}:}{\code{numeric}, convergence code from
      \code{"optim"}.}
    \item{\code{optimOpts}:}{\code{list}; typically just the
      \code{control} arguments passed to \code{\link{optim}}.}
  }
  In addition, the \code{"fitCopula"} class has a slot
  \describe{
    \item{\code{copula}:}{Object of class \code{"\linkS4class{copula}"}.}
  }
  whereas the \code{"fitMvdc"} has
  \describe{
    \item{\code{mvdc}:}{Object of class \code{"mvdc"}.}
  }
}
\references{
  Genest, C., Ghoudi, K., and Rivest, L.-P. (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika} \bold{82}, 543--552.
}
\keyword{classes}
