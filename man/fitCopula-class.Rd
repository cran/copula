\name{fitCopula-class}
\docType{class}
\alias{fitCopula-class}
\alias{fitMvdc-class}
\alias{summaryFitCopula-class}
\alias{summaryFitMvdc-class}
\title{Class "fitCopula" }
\description{ Classes and summaries related to copula model fitting. }
\section{Objects from the Class}{
  Objects can be created by calls to \code{\link{fitCopula}}, \code{\link{fitMvdc}}
  or to their \code{summary} method.
}
\section{Slots}{
  \describe{
    \item{\code{estimate}:}{ \code{numeric}, parameter estimate. }
    \item{\code{var.est}:}{ \code{numeric}, variance matrix estimate of
      the parameter estimator. See note below.  }
    \item{\code{method}:}{ \code{character}, method of estimation. }
    \item{\code{loglik}:}{ \code{numeric}, loglikelihood evaluated at
      the maximizer. }
    \item{\code{convergence}:}{ \code{numeric}, convergence code from
      \code{"optim"}. }
    \item{\code{nsample}:}{ \code{numeric}, integer representing the sample size. }
    \item{\code{copula}:}{ Object of class \code{"\linkS4class{copula}"}. }
    \item{\code{mvdc}:}{ Object of class \code{"mvdc"}. }
  }
}

\references{
  C. Genest, K. Ghoudi and L.-P. Rivest (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika} \bold{82}, 543--552. 
}

\keyword{classes}
