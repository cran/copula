\name{AssocMeasures}
% \docType{methods}
\alias{kendallsTau}
\alias{spearmansRho}
\alias{tailIndex}
\alias{calibKendallsTau}
\alias{calibSpearmansRho}

\alias{kendallsTau-methods}
\alias{kendallsTau,ANY-method}
\alias{kendallsTau,copula-method}
\alias{kendallsTau,archmCopula-method}
\alias{kendallsTau,amhCopula-method}
\alias{kendallsTau,claytonCopula-method}
\alias{kendallsTau,frankCopula-method}
\alias{kendallsTau,gumbelCopula-method}
\alias{kendallsTau,normalCopula-method}
\alias{kendallsTau,tCopula-method}
\alias{kendallsTau,fgmCopula-method}
\alias{kendallsTau,plackettCopula-method}
\alias{kendallsTau,galambosCopula-method}
\alias{kendallsTau,huslerReissCopula-method}
\alias{kendallsTau,tawnCopula-method}
\alias{kendallsTau,tevCopula-method}
\alias{kendallsTau,evCopula-method}

\alias{spearmansRho-methods}
\alias{spearmansRho,ANY-method}
\alias{spearmansRho,copula-method}
\alias{spearmansRho,normalCopula-method}
\alias{spearmansRho,tCopula-method}
\alias{spearmansRho,claytonCopula-method}
\alias{spearmansRho,gumbelCopula-method}
\alias{spearmansRho,frankCopula-method}
\alias{spearmansRho,fgmCopula-method}
\alias{spearmansRho,plackettCopula-method}
\alias{spearmansRho,galambosCopula-method}
\alias{spearmansRho,huslerReissCopula-method}
\alias{spearmansRho,tawnCopula-method}
\alias{spearmansRho,tevCopula-method}
\alias{spearmansRho,evCopula-method}

\alias{tailIndex-methods}
\alias{tailIndex,ANY-method}
\alias{tailIndex,copula-method}
\alias{tailIndex,evCopula-method}
\alias{tailIndex,claytonCopula-method}
\alias{tailIndex,frankCopula-method}
\alias{tailIndex,gumbelCopula-method}
\alias{tailIndex,normalCopula-method}
\alias{tailIndex,tCopula-method}

\alias{calibKendallsTau-methods}
\alias{calibKendallsTau,ANY-method}
\alias{calibKendallsTau,copula-method}
\alias{calibKendallsTau,ellipCopula-method}
\alias{calibKendallsTau,normalCopula-method}
\alias{calibKendallsTau,tCopula-method}
\alias{calibKendallsTau,amhCopula-method}
\alias{calibKendallsTau,claytonCopula-method}
\alias{calibKendallsTau,frankCopula-method}
\alias{calibKendallsTau,gumbelCopula-method}
\alias{calibKendallsTau,plackettCopula-method}
\alias{calibKendallsTau,fgmCopula-method}
\alias{calibKendallsTau,galambosCopula-method}
\alias{calibKendallsTau,huslerReissCopula-method}
\alias{calibKendallsTau,tawnCopula-method}
\alias{calibKendallsTau,tevCopula-method}

\alias{calibSpearmansRho-methods}
\alias{calibSpearmansRho,ANY-method}
\alias{calibSpearmansRho,copula-method}
\alias{calibSpearmansRho,ellipCopula-method}
\alias{calibSpearmansRho,normalCopula-method}
\alias{calibSpearmansRho,tCopula-method}
\alias{calibSpearmansRho,claytonCopula-method}
\alias{calibSpearmansRho,frankCopula-method}
\alias{calibSpearmansRho,gumbelCopula-method}
\alias{calibSpearmansRho,galambosCopula-method}
\alias{calibSpearmansRho,huslerReissCopula-method}
\alias{calibSpearmansRho,tawnCopula-method}
\alias{calibSpearmansRho,tevCopula-method}
\alias{calibSpearmansRho,plackettCopula-method}
\alias{calibSpearmansRho,fgmCopula-method}


\title{Dependence measures for copulas}
\description{
  These functions compute Kendall's Tau, Spearman's Rho, and the tail
  dependence index for bivariate copulas. Calibration is the inverse function: it
  calibrates the copula parameter given the value of Kendall's Tau or Spearman's Rho.
}
\usage{
kendallsTau (copula)
spearmansRho(copula)
tailIndex(copula, \dots)
calibKendallsTau (copula, tau)
calibSpearmansRho(copula, rho)
}
\arguments{
  \item{copula}{an \R object of class \code{"\linkS4class{copula}"}.}
  \item{tau}{ a numerical value of Kendall's Tau in [-1, 1]. }
  \item{rho}{ a numerical value of Spearman's Rho in [-1, 1]. }
  \item{\dots}{ currently nothing. }
}
\details{
  The calibration function in fact returns a moment estimate of the
  parameter for one-parameter copulas.

  When there are no closed-form expressions for Kendall's tau or
  Spearman's rho, the calibration functions use numerical
  approximation techniques (see the last reference). For closed-form
  expressions, see Frees and Valdez (1998). For the t copula, the
  calibration function based on Spearman's rho uses the corresponding
  expression for the normal copula as an approximation.
}
\references{
  E.W. Frees and E.A. Valdez (1998)
  Understanding relationships using copulas.
  \emph{North American Actuarial Journal} \bold{2}, 1--25.

  I. Kojadinovic and J. Yan (2010). Comparison of three semiparametric methods
  for estimating dependence parameters in copula models.
  \emph{Insurance: Mathematics and Economics} \bold{47}, 52--63.
}
\seealso{
  The \code{\linkS4class{acopula}} class objects have slots,
  \code{tau}, \code{lambdaL}, and \code{lambdaU} providing functions for
  \code{kendallsTau()}, and the two tail indices \code{tailIndex()},
  and slot \code{tauInv} for \code{calibKendallsTau()}, see the
  examples, and \code{\link{copGumbel}}, etc.
}
\examples{
gumbel.cop <- gumbelCopula(3)
kendallsTau(gumbel.cop)
spearmansRho(gumbel.cop)
tailIndex(gumbel.cop)

stopifnot(all.equal(kendallsTau(gumbel.cop), copGumbel@tau(3)),

          all.equal(tailIndex(gumbel.cop),
                    c(copGumbel@lambdaL(3), copGumbel@lambdaU(3)),
                    check.attr=FALSE),

          all.equal(calibKendallsTau (gumbel.cop, 0.681),
                    copGumbel@tauInv(0.681))
)


## let us compute the sample versions
x <- rcopula(gumbel.cop, 200)
cor(x, method = "kendall")
cor(x, method = "spearman")
## compare with the true parameter value 3
calibKendallsTau (gumbel.cop, cor(x, method="kendall")[1,2])
calibSpearmansRho(gumbel.cop, cor(x, method="spearman")[1,2])
}
\keyword{multivariate}
