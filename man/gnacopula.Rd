\name{gnacopula}
\alias{gnacopula}
\title{Goodness-of-fit Testing for (Nested) Archimedean Copulas}
\description{
  \code{gnacopula()} conducts a goodness-of-fit test for the given
  (\eqn{H_0}-)copula \code{cop} based on the (copula-)data \code{u}.
}
\usage{
gnacopula(u, cop, n.bootstrap,
          estim.method = eval(formals(enacopula)$method),
          include.K=TRUE, n.MC=0, trafo=c("Hering.Hofert", "Rosenblatt"),
          method=eval(formals(gofTstat)$method), verbose=TRUE, \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of values in \eqn{[0,1]};
    should be (pseudo-/copula-)observations from the copula to be
    tested. Consider applying the function \code{\link{pobs}()} first
    in order to obtain \code{u}.}
  \item{cop}{\eqn{H_0}-\code{"\linkS4class{outer_nacopula}"} with specified
    parameters to be tested for (currently only Archimedean copulas are
    provided).}
  \item{n.bootstrap}{positive integer specifying the number of bootstrap
    replicates.}
  \item{estim.method}{\code{\link{character}}
    string determining the estimation method; see
    \code{\link{enacopula}()}.  We currently only recommend the default
  \code{"mle"} (or maybe \code{"smle"}).}
  \item{include.K}{logical indicating whether the last component, involving the
    Kendall distribution function \code{\link{K}()}, is used in the transformation
    \code{\link{htrafo}()} of Hering and Hofert (2011).  Note that this only
    applies to \code{trafo="Hering.Hofert"}.}
  \item{n.MC}{parameter \code{n.MC} for \code{\link{htrafo}()} (and thus
    for \code{\link{K}()}) if \code{trafo="Hering.Hofert"} and for
    \code{\link{rtrafo}()} if \code{trafo="Rosenblatt"}.}
  \item{trafo}{a \code{\link{character}} string specifying the multivariate
    transformation performed for goodness-of-fit testing, which has to be one
    (or a unique abbreviation) of
    \describe{
      \item{\code{"Hering.Hofert"}}{for the multivariate transformation of
	Hering and Hofert (2011); see \code{\link{htrafo}()}.}
      \item{\code{"Hering.Hofert"}}{for the multivariate transformation of
	Rosenblatt (1952); see \code{\link{rtrafo}()}.}
    }
  }
  \item{method}{a \code{\link{character}} string specifying the
    goodness-of-fit test statistic to be used; see \code{\link{gofTstat}()}.}
  \item{verbose}{if \code{TRUE}, the progress of the bootstrap is
    displayed via \code{\link[utils]{txtProgressBar}}.}
  \item{\dots}{additional arguments passed to \code{\link{enacopula}()}.}
}
\details{
  The function \code{gnacopula()} performs a parametric bootstrap for the
  goodness-of-fit test specified by \code{trafo} and \code{method}.  The
  transformation given by \code{trafo} specifies the multivariate
  transformation which is first applied to the (copula-) data \code{u}
  (typically, the pseudo-observations are used); see
  \code{\link{htrafo}()} or \code{\link{rtrafo}()} for more details. The
  argument \code{method} specifies the particular goodness-of-fit test
  carried out, which is either the Anderson-Darling test for the
  univariate standard uniform distribution (for \code{method="AnChisq"} or
  \code{method="AnGamma"}) in a one-dimensional setup or the tests
  described in Genest et al. (2009) for the
  multivariate standard uniform distribution directly in a multivariate
  setup.  As estimation method, the method provided by
  \code{estim.method} is used.

  Note that a finite-sample correction is made when computing p-values;
  see \code{\link{gofCopula}()} for details.

  A word of warning: Do work carefully with the variety of different
  goodness-of-fit tests that can be performed with \code{gnacopula()}.
  For example, among the possible estimation methods at hand, only
  MLE is known to be consistent (under conditions to be verified).
  Furthermore, for the tests based on the Anderson-Darling test
  statistic, it is theoretically not clear whether the parametric bootstrap
  converges.  Consequently, the results obtained should be treated with
  care.  Moreover, several estimation methods are known to be prone to
  numerical errors (see Hofert et al. (2013)) and are thus not
  recommended to be used in the parametric bootstrap. A warning is given if
  \code{gnacopula()} is called with a method not being MLE.
}
\value{
  \code{gnacopula} returns an \R object of class \code{"htest"}.
  This object contains a list with the bootstrap results including the
  components
  \describe{
    \item{\code{p.value}:}{the bootstrapped p-value;}
    \item{\code{statistic}:}{the value of the test statistic computed
      for the data \code{u};}
    \item{\code{data.name}:}{the name of \code{u};}
    \item{\code{method}:}{a \code{\link{character}} describing the
      goodness-of-fit test applied;}
    \item{\code{estimator}:}{the estimator computed for the data \code{u};}
    \item{\code{bootStats}:}{a list with component \code{estimator}
      containing the estimators for all bootstrap replications and
      component \code{statistic} containing the values of the test statistic
      for each bootstrap replication.}
  }
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

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{gofTstat}()} for the implemented test statistis,
  \code{\link{gtrafo}()} for the multivariate transformation(s)
  \code{\link{htrafo}()} and \code{\link{rtrafo}()} involved and
  \code{\link{K}()} for the Kendall distribution function.

  \code{\link{gofCopula}()} for other (parametric bootstrap) based
  goodness-of-fit tests.
}
\examples{
tau <- 0.5
(theta <- copGumbel@iTau(tau)) # 2
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
uR <- rtrafo(u, cop=copG) # exact
pairs(uR, cex=0.2, gap=.2) # looks good

## plot of the transformed data (Hering and Hofert (2011))
uH <- htrafo(u, cop=copG)
pairs(uH, cex=0.2, gap=.2) # looks good

##' QQ-plot of the supposedly U[0,1] distributed variates
qq1dtrafo <- function(uh, ur) {
    require(lattice)
    d.z <- data.frame(z = c(
                      pchisq(rowSums(qnorm(uh)^2), d), ## Chi^2 -> U[0,1]
                      pgamma(rowSums(-log(uh)), shape=d), ## Gamma -> U[0,1]
                      pchisq(rowSums(qnorm(ur)^2), d),    ## Chi^2 -> U[0,1]
                      pgamma(rowSums(-log(ur)), shape=d)), ## Gamma -> U[0,1]
                      method = rep(c("H.chi2", "H.gam", "R.chi2",
                                     "R.gam"), each=n))
    qqmath(~ z | method, data = d.z, distribution = qunif,
           aspect = "xy", prepanel = prepanel.qqmathline,
           panel = function(x, ...) {
               panel.qqmathline(x, ...)
               panel.qqmath(x, ...)
           })
}
qq1dtrafo(uH, uR)

## what about a meta-Clayton model? (choose larger n.bootstrap in a
## realistic setup)
## note: the parameter of the Clayton copula is only a dummy,
##       it will be estimated anyway
copC <- onacopulaL("Clayton", list(1, 1:d))
res.H0.C <- gnacopula(u, cop=copC, n.bootstrap=10,
	              estim.method="mle")
res.H0.C$p.value # rejection according to 5\% level

## plot of the transformed data (Hering and Hofert (2011)) to see the deviations
## from uniformity
uH <- htrafo(u, cop=copC) # transform the data
pairs(uH, cex=0.2) # clearly visible

## plot of the transformed data (Rosenblatt (1952)) to see the deviations from
## uniformity
uR <- rtrafo(u, cop=copC) # transform the data
pairs(uR, cex=0.2) # clearly visible

## plot the supposedly U[0,1] distributed variates
qq1dtrafo(uH, uR)
## clearly not uniform (as expected)
}
\keyword{htest}
\keyword{distribution}
\keyword{multivariate}
