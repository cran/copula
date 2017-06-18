\name{exchEVTest}
\alias{exchEVTest}
\title{Test of Exchangeability for Certain Bivariate Copulas}
\description{
  Test for assessing the exchangeability of the underlying
  bivariate copula when it is either extreme-value or left-tail
  decreasing.  The test uses the nonparametric estimators of the
  Pickands dependence function studied by Genest and Segers (2009).

  The test statistic is defined in the second reference.
  An approximate p-value for the test statistic is obtained
  by means of a \emph{multiplier} technique if there are no ties in the
  component series of the bivariate data, or by means of an appropriate
  bootstrap otherwise.
}
\usage{
exchEVTest(x, N = 1000,  estimator = c("CFG", "Pickands"),
           ties = NA, ties.method = eval(formals(rank)$ties.method),
           m = 100, derivatives = c("Cn", "An"))
}
\arguments{
  \item{x}{a data matrix that will be transformed to
    pseudo-observations.}
  \item{N}{number of multiplier or boostrap iterations to be used to
    simulate realizations of the test statistic under the null
    hypothesis.}
  \item{estimator}{string specifying which nonparametric estimator of
    the Pickands dependence function \eqn{A()} to use; can be either
    \code{"CFG"} or \code{"Pickands"}; see Genest and Segers (2009).}
  \item{ties}{logical; if \code{FALSE}, approximate p-values are computed by
    means of multiplier bootstrap; if \code{TRUE}, a boostrap adapted to the
    presence of ties in any of the coordinate samples of \code{x}
    is used; the default value of \code{NA} indicates that the
    presence/absence of ties will be checked for automatically.}
   \item{ties.method}{string specifying how ranks should be computed if
     there are ties in any of the coordinate samples of \code{x}; passed
     to \code{\link{pobs}}.}
   \item{derivatives}{a string specifying how the derivatives of the
    unknown copula are estimated; can be either \code{"An"} or \code{"Cn"}.
    The former should be used under the assumption of extreme-value
    dependence.  The latter is faster; see the second reference.}
  \item{m}{integer specifying the size of the integration grid for the
    statistic.}
}
\details{
  More details are available in the references.
}
\value{
  An object of \code{\link{class}} \code{htest} which is a list,
  some of the components of which are
  \item{statistic}{ value of the test statistic. }
  \item{pvalue}{ corresponding approximate p-value. }
}
%\author{Ivan Kojadinovic and Jun Yan}
\references{
  Genest, C. and Segers, J. (2009)
  Rank-based inference for bivariate extreme-value copulas.
  \emph{Annals of Statistics} \bold{37}, 2990--3022.

  Kojadinovic, I. and Yan, J. (2012)
  A nonparametric test of exchangeability for extreme-value and left-tail
  decreasing bivariate copulas.
  \emph{The Scandinavian Journal of Statistics} \bold{39:3}, 480--496.

  Kojadinovic, I. (2017). Some copula inference procedures adapted to
  the presence of ties. \emph{Computational Statistics and Data
    Analysis} \bold{112}, 24--41, \url{http://arxiv.org/abs/1609.05519}.
}
\seealso{ \code{\link{exchTest}}, \code{\link{radSymTest}}, \code{\link{gofCopula}}.
}
\examples{
## Data from an exchangeable left-tail decreasing copulas
exchEVTest(rCopula(200,  gumbelCopula(3)))
\donttest{exchEVTest(rCopula(200, claytonCopula(3)))}

## An asymmetric Khoudraji-Gumbel-Hougaard copula
kc <- khoudrajiCopula(copula1 = indepCopula(),
                      copula2 = gumbelCopula(4),
                      shapes = c(0.4, 0.95))
exchEVTest(rCopula(200, kc))
}
\keyword{htest}
\keyword{multivariate}
