\name{assocMeasures}
% \docType{methods}
\alias{tau}
\alias{rho}
\alias{kendallsTau}
\alias{spearmansRho}
\alias{iTau}
\alias{iRho}
\alias{calibKendallsTau}
\alias{calibSpearmansRho}
\alias{tailIndex}% Deprecated for
\alias{lambda}

\alias{tau-methods}
\alias{tau,ANY-method}
\alias{tau,copula-method}
\alias{tau,archmCopula-method}
\alias{tau,amhCopula-method}
\alias{tau,claytonCopula-method}
\alias{tau,frankCopula-method}
\alias{tau,gumbelCopula-method}
\alias{tau,indepCopula-method}
\alias{tau,joeCopula-method}
\alias{tau,normalCopula-method}
\alias{tau,tCopula-method}
\alias{tau,fgmCopula-method}
\alias{tau,plackettCopula-method}
\alias{tau,galambosCopula-method}
\alias{tau,huslerReissCopula-method}
\alias{tau,tawnCopula-method}
\alias{tau,tevCopula-method}
\alias{tau,evCopula-method}
\alias{tau,acopula-method}
\alias{tau,nacopula-method}
\alias{tau,rotCopula-method}

\alias{rho-methods}
\alias{rho,ANY-method}
\alias{rho,copula-method}
\alias{rho,normalCopula-method}
\alias{rho,tCopula-method}
\alias{rho,claytonCopula-method}
\alias{rho,frankCopula-method}
\alias{rho,gumbelCopula-method}
\alias{rho,amhCopula-method}
\alias{rho,fgmCopula-method}
\alias{rho,plackettCopula-method}
\alias{rho,galambosCopula-method}
\alias{rho,huslerReissCopula-method}
\alias{rho,indepCopula-method}
\alias{rho,tawnCopula-method}
\alias{rho,tevCopula-method}
\alias{rho,evCopula-method}
\alias{rho,acopula-method}
\alias{rho,nacopula-method}
\alias{rho,rotCopula-method}

\alias{lambda-methods}
\alias{lambda,ANY-method}
\alias{lambda,copula-method}
\alias{lambda,evCopula-method}
\alias{lambda,claytonCopula-method}
\alias{lambda,frankCopula-method}
\alias{lambda,gumbelCopula-method}
\alias{lambda,indepCopula-method}
\alias{lambda,amhCopula-method}
\alias{lambda,joeCopula-method}
\alias{lambda,plackettCopula-method}
\alias{lambda,normalCopula-method}
\alias{lambda,tCopula-method}
\alias{lambda,acopula-method}
\alias{lambda,nacopula-method}
\alias{lambda,rotCopula-method}

\alias{iTau-methods}
\alias{iTau,ANY-method}
\alias{iTau,copula-method}
\alias{iTau,ellipCopula-method}
\alias{iTau,normalCopula-method}
\alias{iTau,tCopula-method}
\alias{iTau,amhCopula-method}
\alias{iTau,claytonCopula-method}
\alias{iTau,frankCopula-method}
\alias{iTau,gumbelCopula-method}
\alias{iTau,joeCopula-method}
\alias{iTau,plackettCopula-method}
\alias{iTau,fgmCopula-method}
\alias{iTau,galambosCopula-method}
\alias{iTau,huslerReissCopula-method}
\alias{iTau,tawnCopula-method}
\alias{iTau,tevCopula-method}
\alias{iTau,nacopula-method}
\alias{iTau,acopula-method}
\alias{iTau,rotCopula-method}

\alias{iRho-methods}
\alias{iRho,ANY-method}
\alias{iRho,copula-method}
\alias{iRho,ellipCopula-method}
\alias{iRho,normalCopula-method}
\alias{iRho,tCopula-method}
\alias{iRho,claytonCopula-method}
\alias{iRho,frankCopula-method}
\alias{iRho,gumbelCopula-method}
\alias{iRho,galambosCopula-method}
\alias{iRho,huslerReissCopula-method}
\alias{iRho,tawnCopula-method}
\alias{iRho,tevCopula-method}
\alias{iRho,plackettCopula-method}
\alias{iRho,fgmCopula-method}
\alias{iRho,nacopula-method}
\alias{iRho,rotCopula-method}


\title{Dependence Measures for Bivariate Copulas}
\description{
  These functions compute Kendall's tau, Spearman's rho, and the tail
  dependence index for \emph{bivariate} copulas.  \code{iTau} and \code{iRho},
  sometimes called \dQuote{calibration} functions are the inverses: they
  determine (\dQuote{calibrate}) the copula parameter (which must be
  one-dimensional!) given the value of Kendall's tau or Spearman's rho.
}
\usage{
tau (copula, \dots)
rho (copula, \dots)
lambda(copula, \dots)
iTau (copula, tau, \dots)
iRho (copula, rho, \dots)
}
\arguments{
  \item{copula}{an \R object of class \code{"\linkS4class{copula}"}
    (or also \code{"\linkS4class{acopula}"} or
    \code{"\linkS4class{nacopula}"}; note however that some methods
    may not be available for some copula families).}
  \item{tau}{ a numerical value of Kendall's tau in [-1, 1]. }
  \item{rho}{ a numerical value of Spearman's rho in [-1, 1]. }
  \item{\dots}{ currently nothing. }
}
\details{
  The calibration functions \code{iTau()} and \code{iRho()} in fact
  return a moment estimate of the parameter for one-parameter copulas.

  When there are no closed-form expressions for Kendall's tau or
  Spearman's rho, the calibration functions use numerical
  approximation techniques (see the last reference).  For closed-form
  expressions, see Frees and Valdez (1998).  For the t copula, the
  calibration function based on Spearman's rho uses the corresponding
  expression for the normal copula as an approximation.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  E.W. Frees and E.A. Valdez (1998)
  Understanding relationships using copulas.
  \emph{North American Actuarial Journal} \bold{2}, 1--25.

  Iwan Kojadinovic and Jun Yan (2010)
  Comparison of three semiparametric methods for estimating dependence
  parameters in copula models.
  \emph{Insurance: Mathematics and Economics} \bold{47}, 52--63.
}
\seealso{
  The \code{\linkS4class{acopula}} class objects have slots,
  \code{tau}, \code{lambdaL}, and \code{lambdaU} providing functions for
  \code{tau()}, and the two tail indices \code{lambda()},
  and slot \code{iTau} for \code{iTau()}, see the
  examples and \code{\link{copGumbel}}, etc.
}
\examples{
gumbel.cop <- gumbelCopula(3)
tau(gumbel.cop)
rho(gumbel.cop)
lambda(gumbel.cop)
iTau(joeCopula(), 0.5)
%%--- more "tests" are in ../tests/moments.R
stopifnot(all.equal(tau(gumbel.cop), copGumbel@tau(3)),

          all.equal(lambda(gumbel.cop),
                    c(copGumbel@lambdaL(3), copGumbel@lambdaU(3)),
                    check.attributes=FALSE),

          all.equal(iTau (gumbel.cop, 0.681),
                    copGumbel@iTau(0.681))
)

## let us compute the sample versions
x <- rCopula(200, gumbel.cop)
cor(x, method = "kendall")
cor(x, method = "spearman")
## compare with the true parameter value 3
iTau(gumbel.cop, cor(x, method="kendall" )[1,2])
iRho(gumbel.cop, cor(x, method="spearman")[1,2])
}
\keyword{multivariate}
