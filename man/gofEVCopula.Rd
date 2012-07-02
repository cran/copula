\name{gofEVCopula}
\alias{gofEVCopula}
\title{Goodness-of-fit tests for bivariate extreme-value copulas}
\description{
  Goodness-of-fit tests for extreme-value copulas based on the
  empirical process comparing one of the two nonparameteric
  rank-based estimator of the Pickands dependence function studied in
  Genest and Segers (2009) with a parametric estimate
  of the Pickands dependence function derived under the null hypothesis.
  The test statistic is the Cramer-von Mises functional Sn defined
  in equation (5) of Genest, Kojadinovic, Nešlehová, and Yan (2010).
  Approximate p-values for the test statistic are obtained using a
  parametric bootstrap.
}
\usage{
gofEVCopula(copula, x, N = 1000, method = "mpl",
            estimator = "CFG", m = 1000, verbose = TRUE,
            print.every = NULL, optim.method = "Nelder-Mead")
}
\arguments{
  \item{copula}{ object of class \code{"\linkS4class{evCopula}"}
    representing the hypothesized extreme-value copula family.}
  \item{x}{ a data matrix that will be transformed to pseudo-observations. }
  \item{N}{ number of bootstrap samples to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{method}{ estimation method to be used to estimate the
    dependence parameter(s); can be either \code{"mpl"}
    (maximum pseudo-likelihood), \code{"itau"} (inversion of
    Kendall's tau) or \code{"irho"} (inversion of Spearman's rho). }
  \item{estimator}{ specifies which nonparametric rank-based estimator
    of the unknown Pickands dependence function to use; can be either
    \code{"CFG"} (Caperaa-Fougeres-Genest) or  \code{"Pickands"}.}
  \item{m}{ number of points of the uniform grid on [0,1] used to
    compute the test statistic numerically.}
  \item{print.every}{ is deprecated in favor of \code{verbose}.}
  \item{verbose}{a logical specifying if progress of the bootstrap
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{optim.method}{ the method for \code{"\link{optim}"}. }
}
\details{
  More details can be found in the second reference.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{p.value}{ corresponding approximate p-value. }
  \item{parameter}{ estimates of the parameters for the hypothesized
    copula family. }
}
\note{
  For a given degree of dependence, the most popular extreme-value
  copulas are strikingly similar.
}
\references{
  Genest, C. and Segers, J. (2009). Rank-based inference for bivariate
  extreme-value copulas. \emph{Annals of Statistics} \bold{37}, 2990--3022.

  Genest, C. Kojadinovic, I., Ne\enc{š}{s}lehov\enc{á}{a}, J., and
  Yan, J. (2011). A goodness-of-fit test for bivariate extreme-value copulas.
  \emph{Bernoulli} \bold{17}(1), 253--275.
}
\seealso{ \code{\link{evCopula}}, \code{\link{evTestC}}, \code{\link{evTestA}},
  \code{\link{evTestK}}, \code{\link{gofCopula}}, \code{\link{An}}.
}
\examples{
\dontrun{% because the bootstraps are slow, particularly for the non-Gumbel ones
x <- rCopula(100, claytonCopula(3))

## Does the Gumbel family seem to be a good choice?
gofEVCopula(gumbelCopula(1), x)

## The same with different estimation methods
gofEVCopula(gumbelCopula(1), x, method="itau")
gofEVCopula(gumbelCopula(1), x, method="irho")

## The same with different extreme-value copulas
gofEVCopula(galambosCopula(1), x)
gofEVCopula(galambosCopula(1), x, method="itau")
gofEVCopula(galambosCopula(1), x, method="irho")

gofEVCopula(huslerReissCopula(1), x)
gofEVCopula(huslerReissCopula(1), x, method="itau")
gofEVCopula(huslerReissCopula(1), x, method="irho")

gofEVCopula(tevCopula(0, df.fixed=TRUE), x)
gofEVCopula(tevCopula(0, df.fixed=TRUE), x, method="itau")
gofEVCopula(tevCopula(0, df.fixed=TRUE), x, method="irho")
}
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
