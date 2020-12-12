\newcommand{\CRANpkg}{\href{https://CRAN.R-project.org/package=#1}{\pkg{#1}}}
%%
\name{copula-package}
\alias{copula-package}
\docType{package}
\title{Multivariate Dependence Modeling with Copulas}
\description{
  The \pkg{copula} package provides (S4) classes of commonly used
  elliptical, (nested) Archimedean, extreme value and other copula families;
  methods for density, distribution, random number generation, and plots.

  Fitting copula models and goodness-of-fit tests.
  Independence and serial (univariate and multivariate) independence tests,
  and other copula related tests.
}
\details{

  The DESCRIPTION file:
  \packageDESCRIPTION{copula}
  \packageIndices{copula}

  The \pkg{copula} package provides
  \itemize{
    \item Classes (S4) of commonly used copulas including
    elliptical (normal and t; \code{\linkS4class{ellipCopula}}),
    Archimedean (Clayton, Gumbel, Frank, Joe, and Ali-Mikhail-Haq; ;
    \code{\linkS4class{archmCopula}} and \code{\linkS4class{acopula}}),
    extreme value (Gumbel, Husler-Reiss, Galambos, Tawn, and t-EV; \code{\linkS4class{evCopula}}),
    and other families (Plackett and Farlie-Gumbel-Morgenstern).
    \item Methods for density, distribution, random number generation
    (\code{\link{dCopula}}, \code{\link{pCopula}} and \code{\link{rCopula}});
    bivariate dependence measures (\code{\link{rho}}, \code{\link{tau}},
    etc), perspective and contour plots.
    \item Functions (and methods) for fitting copula models including
    variance estimates (\code{\link{fitCopula}}).
%%% FIXME: Add more \code{\link{...}} entries:
    \item Independence tests among random variables and vectors.
    \item Serial independence tests for univariate and multivariate
    continuous time series.
    \item Goodness-of-fit tests for copulas based on multipliers, and
    the parametric bootstrap, with several transformation options.
    \item Bivariate and multivariate tests of extreme-value dependence.
    \item Bivariate tests of exchangeability.
  }%% FIXME 2: Eventually merge the two itemized lists
  Now with former package \pkg{nacopula} for working with nested Archimedean copulas.
  Specifically,
  \itemize{
    \item it provides procedures for computing function values and cube
    volumes (\code{\link{prob}}),
    \item characteristics such as Kendall's tau and tail dependence
    coefficients (via family objects, e.g.,
    \code{\link[=acopula-families]{copGumbel}}),
    \item efficient sampling algorithms (\code{\link{rnacopula}}),
    \item various estimators and goodness-of-fit tests.
    \item The package also contains related univariate distributions and special functions
    such as the Sibuya distribution (\code{\link{Sibuya}}), the
    polylogarithm (\code{\link{polylog}}), Stirling and Eulerian numbers
    (\code{\link{Eulerian}}).
  }

  Further information is available in the following \link{vignettes}:
  \tabular{ll}{
    \code{nacopula-pkg} \tab Nested Archimedean Copulas Meet R (\url{../doc/nacopula-pkg.pdf})\cr
    \code{Frank-Rmpfr}  \tab {Numerically Stable Frank via Multiprecision in R} (\url{../doc/Frank-Rmpfr})\cr
  }

  For a list of exported functions, use \code{help(package = "copula")}.
}
%% \author{
%%   Marius Hofert, % \email{marius.hofert@math.ethz.ch}
%%   Ivan Kojadinovic, % \email{ivan.kojadinovic@univ-pau.fr}
%%   Martin Maechler, and % \email{maechler@stat.math.ethz.ch},
%%   Jun Yan.% \email{jun.yan@uconn.edu}.

%%   Maintainer: Currently,
%%   Martin Maechler \email{maechler@stat.math.ethz.ch}.
%% }
\references{
  Yan, J. (2007)
  Enjoy the Joy of Copulas: With a Package \pkg{copula}.
  \emph{Journal of Statistical Software} \bold{21}(4), 1--21.
  \url{https://www.jstatsoft.org/v21/i04/}.

  Kojadinovic, I. and Yan, J. (2010).
  Modeling Multivariate Distributions with Continuous Margins Using the
  copula R Package.
  \emph{Journal of Statistical Software} \bold{34}(9), 1--20.
  \url{https://www.jstatsoft.org/v34/i09/}.

  Hofert, M. and \enc{Mächler}{Maechler}, M. (2011),
  Nested Archimedean Copulas Meet R: The nacopula Package.,
  \emph{Journal of Statistical Software} \bold{39}(9), 1--20.
  \url{https://www.jstatsoft.org/v39/i09/}.

  Nelsen, R. B. (2006)
  \emph{An introduction to Copulas}. Springer, New York.
}
\seealso{
  The following CRAN packages currently use (\sQuote{depend on}) \pkg{copula}:
  \CRANpkg{CoClust}, \CRANpkg{copulaedas}, \CRANpkg{Depela},
  \CRANpkg{HAC}, \CRANpkg{ipptoolbox}, \CRANpkg{vines}.
  %% and \CRANpkg{nacopula}, but that's only formal
}
\examples{
## Some of the more important functions (and their examples) are
\donttest{
example(fitCopula)## fitting Copulas
example(fitMvdc)  ## fitting multivariate distributions via Copulas
example(nacopula) ## nested Archimedean Copulas

## Independence Tests:  These also draw a 'Dependogram':
example(indepTest)       ## Testing for Independence
example(serialIndepTest) ## Testing for Serial Independence
}% not in R CMD check
}
\keyword{package}
