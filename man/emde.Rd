\name{emde}
\alias{emde}
\title{Minimum Distance Estimators for (Nested) Archimedean Copulas}
\description{
  Compute minimum distance estimators for (nested) Archimedean copulas.
}
\usage{
emde(u, cop,
     method=c("mde.chisq.CvM", "mde.chisq.KS",
              "mde.gamma.CvM", "mde.gamma.KS"),
     interval=initOpt(cop@copula@name),
     include.K = FALSE, repara = TRUE, \dots)
}
\arguments{
  \item{u}{\eqn{n\times d}{n x d}-matrix of (pseudo-)observations (each
    value in \eqn{[0,1]}) from the copula, where \eqn{n} denotes the
    sample size and \eqn{d} the dimension.}
  \item{cop}{\code{\linkS4class{outer_nacopula}} to be estimated
    (currently only Archimedean copulas are provided).}
  \item{method}{a \code{\link{character}} string specifying the distance
    method, which has to be one (or a unique abbreviation) of
    \describe{
      \item{\code{"mde.chisq.CvM"}}{map to an Erlang distribution and
	using a chi-square distribution and Cram\enc{é}{e}r-von Mises distance;}
      \item{\code{"mde.chisq.KS"}}{map to an Erlang distribution and
	using a chi-square distribution and Kolmogorov-Smirnov distance;}
      \item{\code{"mde.gamma.CvM"}}{map to an Erlang distribution and
	using a Erlang distribution and Cram\enc{é}{e}r-von Mises distance;}
      \item{\code{"mde.gamma.KS"}}{map to an Erlang distribution and using
	a \ifelse{latex}{Kolmogorov-Smir-nov}{Kolmogorov-Smirnov} distance.}

    }
    The four methods are described in Hofert et al. (2011a); see also the
    \sQuote{Details} section.}
  \item{interval}{bivariate vector denoting the interval where
    optimization takes place.  The default is computed as described in Hofert et al.(2011a).}
  \item{include.K}{logical indicating whether the last component, the
    (possibly numerically challenging) Kendall distribution function \eqn{K}, is used
    (\code{include.K=TRUE}) or not.  Note that the default is
    \code{\link{FALSE}} here, where it is \code{TRUE} in the underlying
    \code{\link{htrafo}()} function.}
  \item{repara}{logical indicating whether the distance function to be
    optimized is reparametrized (the default); see the code for more
    details.}
  \item{\dots}{additional arguments passed to \code{\link{optimize}}.}
}
\details{
   First, \code{\link{htrafo}} is applied to map the \eqn{n\times d}{n x d}-matrix
   of given realizations to a \eqn{n\times d}{n x d}-matrix or
   \eqn{n\times (d-1)}{n x (d-1)}-matrix, depending on whether the last
   component is included (\code{include.K=TRUE}) or not.  Second, using
   either the sum of squares of the standard normal quantile function
   (\code{method="mde.chisq.CvM"} and
   \code{method="mde.chisq.KS"}) or the sum of negative logarithms
   (\code{method="mde.gamma.CvM"} and \code{method="mde.gamma.KS"}), a map to a
   chi-square or an Erlang distribution is applied, respectively.  Finally, a
   Cram\enc{é}{e}r-von Mises (\code{method="mde.chisq.CvM"} and
   \code{method="mde.gamma.CvM"}) or Kolmogorov-Smirnov
   (\code{method="mde.chisq.KS"} and \code{method="mde.gamma.KS"}) distance is
   applied.  This is repeated in an optimization until the copula parameter is
   found such that this distance is minimized.

   Note that the same transformations as described above are applied for
   goodness-of-fit testing; see the \sQuote{See Also} section).
}
\value{
  \code{\link{list}} as returned by \code{\link{optimize}}, including the
  minimum distance estimator.
}
\author{Marius Hofert}
\references{
   Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2011a),
   Estimators for Archimedean copulas in high dimensions: A comparison,
   to be submitted.

   Hering, C. and Hofert, M. (2011),
	Goodness-of-fit tests for Archimedean copulas in large dimensions,
	submitted.
}
\seealso{
  \code{\link{enacopula}} (wrapper for different estimators),
  \code{\link{gnacopula}} (wrapper for different goodness-of-fit tests),
  \code{\link{htrafo}} (transformation to a multivariate uniform
   distribution), and
  \code{\link{K}} (Kendall distribution function).
}
\examples{
tau <- 0.25
(theta <- copGumbel@tauInv(tau)) # 4/3
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

set.seed(1)
n <- 200
U <- rnacopula(n, cop)

(meths <- eval(formals(emde)$method)) # "mde.chisq.CvM", ...
fun <- function(meth, u, cop, theta){
	run.time <- system.time(val <- emde(u, cop=cop, method=meth)$minimum)
	list(value=val, error=val-theta, utime.ms=1000*run.time[[1]])
}
(res <- sapply(meths, fun, u=U, cop=cop, theta=theta))
}
\keyword{models}
