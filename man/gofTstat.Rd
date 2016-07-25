\name{gofTstat}
\alias{gofTstat}
\title{Goodness-of-fit Test Statistics}
\description{
  \code{gofTstat()} computes various goodness-of-fit test statistics
  typically used in \code{\link{gofCopula}(*, simulation = "pb")}.
}
\usage{
gofTstat(u, method = c("Sn", "SnB", "SnC", "AnChisq", "AnGamma"),
         useR = FALSE, \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of values in \eqn{[0,1]},
    supposedly independent uniform observations in the
    hypercube, that is, \eqn{U_i \sim U[0,1]^d}{U_i ~ U[0,1]^d}, i.i.d.,
    for \eqn{i \in \{1,\dots,n\}}{i in 1..n}.}
  \item{method}{a \code{\link{character}} string specifying the
    goodness-of-fit test statistic to be used, which has to be one (or a unique
    abbreviation) of
    \describe{
      \item{\code{"Sn"}}{for computing the test statistic
	\eqn{S_n} from Genest,
	R\enc{é}{e}millard, Beaudoin (2009).}
      \item{\code{"SnB"}}{for computing the test statistic
	\eqn{S_n^{(B)}}{S_n^(B)} from Genest,
	R\enc{é}{e}millard, Beaudoin (2009).}
      \item{\code{"SnC"}}{for computing the test statistic
	\eqn{S_n^{(C)}}{S_n^(C)} from Genest et al. (2009).}
      \item{\code{"AnChisq"}}{Anderson-Darling test statistic for
	computing (supposedly) \eqn{\mathrm{U}[0,1]}{U[0,1]}-distributed
	(under \eqn{H_0}) random variates via the distribution function
	of the chi-square distribution with \eqn{d}
	degrees of freedom. To be more precise, the Anderson-Darling
	test statistc of the variates
	\deqn{\chi_d^2\Bigl(\sum_{j=1}^d(\Phi^{-1}(u_{ij}))^2\Bigr)}{%
	  pchisq((Phi^{-1}(u_{i1}))^2+...+(Phi^{-1}(u_{id}))^2, df=d)}
	is computed (via \code{ADGofTest::ad.test}), where
	\eqn{\Phi^{-1}}{Phi^{-1}} denotes the quantile function
	of the standard normal distribution function,
	\eqn{\chi_d^2}{pchisq(.,df=d)} denotes the distribution function
	of the chi-square distribution with \eqn{d} degrees of freedom,
	and \eqn{u_{ij}} is the \eqn{j}th component in the \eqn{i}th row
	of \code{u}.}
      \item{\code{"AnGamma"}}{similar to \code{method="AnChisq"} but
	based on the variates
	\deqn{\Gamma_d\Bigl(\sum_{j=1}^d(-\log u_{ij})\Bigr),}{%
	  pgamma(-log(u_{i1})-...-log(u_{id}), shape=d),}
	where \eqn{\Gamma_d}{pgamma(.,shape=d)} denotes the distribution
	function of the gamma distribution with shape parameter \eqn{d}
	and shape parameter one	(being equal to an
	Erlang(\eqn{d}) distribution function).}
  }}
  \item{useR}{logical indicating whether an \R or C implementation is used.}
  \item{\dots}{additional arguments passed for computing the different
    test statistics.}
}
\details{
  This function should be used with care.  The different test statistics
  were implemented (partly) for different purposes and goodness-of-fit
  tests and should be used only with knowledge about such (see the
  references for more details).
}
\value{
  The value of the test statistic, a \code{\link{numeric}}.
}
%\author{Marius Hofert and Martin Maechler}
\references{
  Genest, C., R\enc{é}{e}millard, B., and Beaudoin, D. (2009),
  Goodness-of-fit tests for copulas: A review and a power study
  \emph{Insurance: Mathematics and Economics} \bold{44}, 199--213.

  Rosenblatt, M. (1952),
  Remarks on a Multivariate Transformation,
  \emph{The Annals of Mathematical Statistics} \bold{23}, 3, 470--472.

  Hering, C. and Hofert, M. (2014),
  Goodness-of-fit tests for Archimedean copulas in high dimensions,
  \emph{Innovations in Quantitative Risk Management}.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.
}
\seealso{
  \code{\link{gofCopula}()} for goodness-of-fit tests where (some of)
  these test statistics are used.
}
\examples{
## generate data
cop <- archmCopula("Gumbel", param=iTau(gumbelCopula(), 0.5), dim=5)
set.seed(1)
U <- rCopula(1000, cop)

## compute Sn (as is done in a parametric bootstrap, for example)
Uhat <- pobs(U) # pseudo-observations
u <- cCopula(Uhat, copula = cop) # Rosenblatt transformed data (with correct copula)
gofTstat(u, method = "Sn", copula = cop) # compute test statistic Sn; requires copula argument
}
\keyword{htest}
\keyword{goodness-of-fit}
\keyword{distribution}
\keyword{multivariate}
