\name{generator}
% \docType{methods}
\alias{genFun}
\alias{genFun-methods}
\alias{genFun,claytonCopula-method}
\alias{genFun,frankCopula-method}
\alias{genFun,gumbelCopula-method}
\alias{genFun,amhCopula-method}

\alias{genInv}
\alias{genInv-methods}
\alias{genInv,claytonCopula-method}
\alias{genInv,frankCopula-method}
\alias{genInv,gumbelCopula-method}
\alias{genInv,amhCopula-method}

\alias{genFunDer1}
\alias{genFunDer1-methods}
\alias{genFunDer1,claytonCopula-method}
\alias{genFunDer1,frankCopula-method}
\alias{genFunDer1,gumbelCopula-method}
\alias{genFunDer1,amhCopula-method}

\alias{genFunDer2}
\alias{genFunDer2-methods}
\alias{genFunDer2,claytonCopula-method}
\alias{genFunDer2,frankCopula-method}
\alias{genFunDer2,gumbelCopula-method}
\alias{genFunDer2,amhCopula-method}

\alias{Afun}
\alias{Afun-methods}
\alias{Afun,galambosCopula-method}
\alias{Afun,gumbelCopula-method}
\alias{Afun,huslerReissCopula-method}

\title{Generator functions for Archimedean and extreme value copulas}
\description{
 Methods to evaluate the generator function, the inverse generator
 function, and derivatives of the generator function.
}
\usage{
genFun(copula, u)
genInv(copula, s)
genFunDer1(copula, u)
genFunDer2(copula, u)
Afun(copula, w)
}
\arguments{
  \item{copula}{an object of class \code{"copula"}. }
  \item{u, s, w}{numerical vector at which these functions are to be
    evaluated. }
}
\details{
  \code{"genFun"} and \code{"genInv"} are, respectively, the generator
  function and its inverse function for an Archimedean copula.

  \code{"genFunDer1"} and \code{"genFunDer2"} are, respectively, the 1st
  and 2nd derivatives of the generator function for an Archimedean copula.

  \code{"Afun"} is the generator function of an extreme value copula.
}
\keyword{methods}
\keyword{multivariate}
