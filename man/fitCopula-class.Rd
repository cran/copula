\name{fitCopula-class}
\docType{class}
\alias{fitCopula-class}
\alias{fitMvdc-class}
\alias{summaryFitCopula-class}
\alias{summaryFitMvdc-class}
\title{Class "fitCopula" }
\description{ Classes and summaries related to copula model fitting. }
\section{Objects from the Class}{
  Objects can be created by calls to \code{fitCopula}, \code{fitMvdc}
  or to their \code{summary} method.
}
\section{Slots}{
  \describe{
    \item{\code{est}:}{ \code{numeric}, parameter estimate. }
    \item{\code{var.est}:}{ \code{numeric}, variance matrix estimate of
      the parameter estimator. See note below.  }
    \item{\code{loglik}:}{ \code{numeric}, loglikelihood evaluated at
      the maximizer. }
    \item{\code{convergence}:}{ \code{numeric}, convergence code from
      \code{"optim"}. }
    \item{\code{nsample}:}{ \code{numeric}, integer representing the sample size. }
    \item{\code{copula}:}{ Object of class \code{"copula"}. }
    \item{\code{mvdc}:}{ Object of class \code{"mvdc"}. }
  }
}

\note{
  The variance matrix of the estimate, \code{"var.est"}, in the returned
  object is only valid for full likelihood maximization. When
  pseudo-observations are used, they under-estimate the variation. See
  Genest, Ghoudi, and Rivest (1995).
}

\references{
  C. Genest, K. Ghoudi and L.-P. Rivest (1995). A semiparametric
  estimation procedure of dependence parameters in multivariate
  families of distributions. \emph{Biometrika}, 82, 543-552. 

  J. Yan (2006), Multivariate Modeling with Copulas and Engineering
  Applications. In \emph{Handbook of Engineering Statistics}, Ed. Pham,
  Springer.
}

\keyword{classes}
