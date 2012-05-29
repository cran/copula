\name{gnacopula}
\alias{gnacopula}
\alias{gtrafouni}
\title{Goodness-of-Fit Testing for (Nested) Archimedean Copulas}
\description{
  \code{gnacopula()} conducts a goodness-of-fit test for the given
  (\eqn{H_0}-)copula \code{cop} based on the (copula-)data \code{u}.

  \code{gtrafouni()} transforms supposedly \eqn{\mathrm{U}[0,1]^d}{U[0,1]^d}
  distributed vectors of random variates to univariate data (for testing in a
  one-dimensional setup).
}
\usage{
gnacopula(u, cop, n.bootstrap,
          estimation.method = eval(formals(enacopula)$method),
          include.K=TRUE, n.MC=0, trafo= c("Hering.Hofert", "Rosenblatt"),
          method=eval(formals(gtrafouni)$method), verbose=TRUE, \dots)
gtrafouni(u, method = c("chisq", "gamma", "Remillard", "Genest"))
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-/copula-)observations (each
    value in \eqn{[0,1]}) from the copula to be tested.  Consider applying the
    function \code{\link{pobs}} first in order to obtain \code{u}.}
  \item{cop}{\eqn{H_0}-\code{"\linkS4class{outer_nacopula}"} with specified
    parameters to be tested for (currently only Archimedean copulas are
    provided).}
  \item{n.bootstrap}{positive integer specifying the number of bootstrap replicates.}
  \item{estimation.method}{\code{\link{character}}
    string determining the estimation method; see
    \code{\link{enacopula}}.  We currently only recommend the default
  \code{"mle"} (or maybe \code{"smle"}).}
  \item{include.K}{logical indicating whether the last component, involving the
    Kendall distribution function \code{\link{K}}, is used in the transformation
    \code{\link{htrafo}} of Hering and Hofert (2011).  Note that this only
    applies to \code{trafo="Hering.Hofert"}.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{htrafo}} (and thus
    for \code{\link{K}}) if \code{trafo="Hering.Hofert"} and for
    \code{\link{rtrafo}} if \code{trafo="Rosenblatt"}.}
  \item{trafo}{a \code{\link{character}} string specifying the multivariate
    transformation performed for goodness-of-fit testing, which has to be one
    (or a unique abbreviation) of
    \describe{
      \item{\code{"Hering.Hofert"}}{for the multivariate transformation of
	Hering and Hofert (2011); see \code{\link{htrafo}}.}
      \item{\code{"Hering.Hofert"}}{for the multivariate transformation of
	Rosenblatt (1952); see \code{\link{rtrafo}}.}
    }
  }
  \item{method}{a \code{\link{character}} string specifying the
    goodness-of-fit method to be used, which has to be one (or a unique
    abbreviation) of
    \describe{
      \item{\code{"chisq"}}{for computing (supposedly)
	\eqn{\mathrm{U}[0,1]}{U[0,1]}-distributed (under \eqn{H_0}) random variates
	via the distribution function of a chi-square distribution with \eqn{d}
	degrees of freedom. To be more precise, the variates
	\deqn{\chi_d^2\Bigl(\sum_{j=1}^d(\Phi^{-1}(u_{ij}))^2\Bigr)}{%
	  pchisq((Phi^{-1}(u_{i1}))^2+...+(Phi^{-1}(u_{id}))^2, df=d)}
	are returned, where \eqn{\Phi^{-1}}{Phi^{-1}} denotes the quantile function
	of the standard normal distribution function, \eqn{\chi_d^2}{pchisq(.,df=d)}
	denotes the distribution function of the chi-square distribution with \eqn{d}
	degrees of freedom, and \eqn{u_{ij}} is the \eqn{j}th component
	in the \eqn{i}th row of \code{u}.}

      \item{\code{"gamma"}}{similar to \code{method="chisq"} but computing
	\deqn{\Gamma_d\Bigl(\sum_{j=1}^d(-\log u_{ij})\Bigr),}{%
	  pgamma(-log(u_{i1})-...-log(u_{id}), shape=d),}
	where \eqn{\Gamma_d}{pgamma(.,shape=d)} denotes the distribution function of
	the gamma distribution with shape parameter \eqn{d} and shape parameter one
	(being equal to an Erlang(\eqn{d}) distribution function).}
      \item{\code{"Remillard"}}{for computing the test statistic
	\eqn{\mathrm{S_n^{(B)}}}{S[n]^(B)} from Genest, R\enc{é}{e}millard, Beaudoin (2009).}
      \item{\code{"Genest"}}{for computing the test statistic
	\eqn{\mathrm{S_n^{(C)}}}{S[n]^(C)} from Genest et al.(2009).}
  }}
  \item{verbose}{if \code{TRUE}, the progress of the bootstrap is
    displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{\dots}{additional arguments to \code{\link{enacopula}}.}
}
\details{
  The function \code{gnacopula} performs a bootstrap for the
  goodness-of-fit test specified by \code{trafo} and \code{method}.  The
  transformation given by \code{trafo} specifies the multivariate
  transformation which is first applied to the (copula-) data \code{u}
  (typically, the pseudo-observations are used); see
  \code{\link{htrafo}} or \code{\link{rtrafo}} for more details. The
  argument \code{method} specifies the particular goodness-of-fit test
  carried out, which is either the Anderson-Darling test for the
  univariate standard uniform distribution (for \code{method="chisq"} or
  \code{method="gamma"}) in a one-dimensional setup or the tests
  described in Genest, R\enc{é}{e}millard, Beaudoin (2009) for the
  multivariate standard uniform distribution directly in a multivariate
  setup.  As estimation method, the method provided by
  \code{estimation.method} is used.

  A word of warning: Do work carefully with the variety of different
  goodness-of-fit tests that can be performed with \code{gnacopula()}.
  For example, among the possible estimation methods at hand, only
  MLE is known to be consistent (under conditions to be verified).
  Furthermore, for the tests based on the Anderson-Darling test
  statistic, it is theoretically not clear whether the bootstrap
  converges.  Consequently, the results obtained should be treated with
  care.  Moreover, several estimation methods are known to be prone to
  numerical errors (see Hofert et al. (2011a)) and are thus not
  recommended to be used in the bootstrap. A warning is given if
  \code{gnacopula()} is called with a method not being MLE.

  Since \R is widely used by practitioners, a word of warning concerning
  goodness-of-fit tests \emph{in general} is also advisable.
  Goodness-of-fit tests are often (ab)used in practice to
  \dQuote{justify} an assumption under which one then continues to work
  (carelessly).  From a mathematical point of view, this is not correct.
}
\value{
  \code{gnacopula} returns an \R object of class \code{"htest"}.
  This object contains a list with the bootstrap results including the components
  \describe{
    \item{\code{p.value}:}{the bootstrapped p-value;}
    \item{\code{statistic}:}{the value of the test statistic computed for the data \code{u};}
    \item{\code{data.name}:}{the name of \code{u};}
    \item{\code{method}:}{a \code{\link{character}} describing the
      goodness-of-fit test applied;}
    \item{\code{estimator}:}{the estimator computed for the data \code{u};}
    \item{\code{bootStats}:}{a list with component \code{estimator}
      containing the estimators for all bootstrap replications and
      component \code{statistic} containing the values of the test statistic
      for each bootstrap replication.}
  }

  \code{gtrafouni} returns a numeric vector of length \eqn{n}
  (=\code{nrow(u)}) containing the univariate transformed values.
}
\author{Marius Hofert, Martin Maechler.}
\references{
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009),
  Goodness-of-fit tests for copulas: A review and a power study
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--213.

  Rosenblatt, M. (1952),
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics} \bold{23}, 3, 470--472.

  Hering, C. and Hofert, M. (2011),
  Goodness-of-fit tests for Archimedean copulas in large dimensions,
  submitted.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011b),
  Likelihood inference for Archimedean copulas,
  submitted.
}
\seealso{
  \code{\link{gtrafo}} for the multivariate transformation(s)
  \code{\link{htrafo}} and \code{\link{rtrafo}} involved and
  \code{\link{K}} for the Kendall distribution function.

  \code{\link{gofCopula}} for other (parametric bootstrap) based
  GOF tests.
}
\examples{
tau <- 0.5
(theta <- copGumbel@tauInv(tau)) # 2
d <- 5
(copG <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 1000
x <- rnacopula(n, copG)
x <- qnorm(x) # x now follows a meta-Gumbel model with N(0,1) marginals
u <- pobs(x) # build pseudo-observations

## check if the data comes from a meta-Gumbel model (choose larger n.bootstrap
## in a realistic setup)
res.H0.G <- gnacopula(u, cop=copG, n.bootstrap=10)
## => uses the transformation of Hering and Hofert (2011), including
##    the Kendall distribution function K and the mapping to a univariate
##    setting via the chi-square distribution. The final test carried out
##    is the Anderson-Darling test.
res.H0.G$p.value # non-rejection according to 5\% level

## plot of the transformed data (Rosenblatt (1952))
u.prime <- rtrafo(u, cop=copG) # exact
pairs(u.prime, cex=0.2) # looks good

## plot of the transformed data (Hering and Hofert (2011))
u.prime. <- htrafo(u, cop=copG)
pairs(u.prime., cex=0.2) # looks good

## what about a meta-Clayton model? (choose larger n.bootstrap in a
## realistic setup)
## note: the parameter of the Clayton copula is only a dummy,
##       it will be estimated anyway
copC <- onacopulaL("Clayton", list(1, 1:d))
res.H0.C <- gnacopula(u, cop=copC, n.bootstrap=10,
	              estimation.method="mle")
res.H0.C$p.value # rejection according to 5\% level

## plot of the transformed data (Hering and Hofert (2011)) to see the deviations
## from uniformity
u.prime <- htrafo(u, cop=copC) # transform the data
pairs(u.prime, cex=0.2) # clearly visible

## plot of the transformed data (Rosenblatt (1952)) to see the deviations from
## uniformity
u.prime. <- rtrafo(u, cop=copC) # transform the data
pairs(u.prime., cex=0.2) # clearly visible

## plot the supposedly U[0,1] distributed variates
z <- gtrafouni(u.prime)
plot(1:length(z), z)
## a bit harder to see, but not perfectly uniform (as expected)
}
\keyword{htest}
\keyword{distribution}
\keyword{multivariate}
