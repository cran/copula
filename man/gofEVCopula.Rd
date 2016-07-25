\name{gofEVCopula}
\alias{gofEVCopula}
\title{Goodness-of-fit Tests for Bivariate Extreme-Value Copulas}
\description{
  Goodness-of-fit tests for extreme-value copulas based on the
  empirical process comparing one of the two nonparameteric
  rank-based estimator of the Pickands dependence function studied in
  Genest and Segers (2009) with a parametric estimate
  of the Pickands dependence function derived under the null hypothesis.
  The test statistic is the Cramer-von Mises functional Sn defined
  in Equation (5) of Genest, Kojadinovic, G. \enc{Nešlehová}{Neslehova}, and Yan (2010).
  Approximate p-values for the test statistic are obtained using a
  parametric bootstrap.
}
\usage{
gofEVCopula(copula, x, N = 1000,
            method = c("mpl", "ml", "itau", "irho"),
            estimator = c("CFG", "Pickands"), m = 1000,
            verbose = interactive(),
            optim.method = "BFGS")
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
  \item{verbose}{a logical specifying if progress of the bootstrap
    should be displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{optim.method}{a string, the \code{method} for
    \code{\link{optim}()}.  In \pkg{copula} versions \code{0.999-14} and
    earlier, the default was \code{"Nelder-Mead"}, but now is set to the
    same as for \code{\link{fitCopula}()}.}
%% Notably, that gave tons of warnings
%%  one-dimensional optimization by Nelder-Mead is unreliable: use "Brent" or optimize() directly
}
\details{
  More details can be found in the second reference.

  The former argument \code{print.every} is deprecated and not
  supported anymore; use \code{verbose} instead.
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
  For a given degree of dependence, the most popular bivariate
  extreme-value copulas are strikingly similar.
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Genest, C. and Segers, J. (2009). Rank-based inference for bivariate
  extreme-value copulas. \emph{Annals of Statistics} \bold{37}, 2990--3022.

  Genest, C. Kojadinovic, I., G. \enc{Nešlehová}{Neslehova}, J., and
  Yan, J. (2011). A goodness-of-fit test for bivariate extreme-value copulas.
  \emph{Bernoulli} \bold{17}(1), 253--275.
}
\seealso{ \code{\link{evCopula}}, \code{\link{evTestC}}, \code{\link{evTestA}},
  \code{\link{evTestK}}, \code{\link{gofCopula}}, \code{\link{An}}.
}
\examples{
n <- 100; N <- 1000 # realistic (but too large currently for CRAN checks)
n <-  60; N <-  200 # (for now)
x <- rCopula(n, claytonCopula(3))

\donttest{% because the bootstraps are slow,..
## Does the Gumbel family seem to be a good choice?
gofEVCopula(gumbelCopula(), x)
}

## The same with different (and cheaper) estimation methods:
gofEVCopula(gumbelCopula(), x, N=N, method="itau")
gofEVCopula(gumbelCopula(), x, N=N, method="irho")

\donttest{% non-gumbel bootstraps are slow
## The same with different extreme-value copulas
gofEVCopula(galambosCopula(), x)
gofEVCopula(galambosCopula(), x, method="itau")
gofEVCopula(galambosCopula(), x, method="irho")

gofEVCopula(huslerReissCopula(), x)
gofEVCopula(huslerReissCopula(), x, method="itau")
gofEVCopula(huslerReissCopula(), x, method="irho")

gofEVCopula(tevCopula(df.fixed=TRUE), x)
gofEVCopula(tevCopula(df.fixed=TRUE), x, method="itau")
gofEVCopula(tevCopula(df.fixed=TRUE), x, method="irho")
}% dont..
}
\keyword{htest}
\keyword{models}
\keyword{multivariate}
