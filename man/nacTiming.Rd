\name{nacFrail.time}
\alias{nacFrail.time}
\title{Timing for Sampling Frailties of Nested Archimedean Copulas}
\usage{
nacFrail.time(n, family, taus, digits = 3, verbose = FALSE)
}
\description{
  This function provides measurements of user run times for the frailty
  variables involved in a nested Archimedean copula.
}
\arguments{
  \item{n}{integer specifying the sample size to be used for the random
    variates \eqn{V_0}{V0} and \eqn{V_{01}}{V01}.}
  \item{family}{the Archimedean family (class
    \code{"\linkS4class{acopula}"}) for which \eqn{V_0}{V0} and
    \eqn{V_{01}}{V01} are sampled.}
  \item{taus}{\code{\link{numeric}} vector of Kendall's taus.  This vector is
    converted to a vector of copula parameters \eqn{theta}, which then serve as
    \eqn{theta_0}{theta0} and \eqn{theta_1}{theta1} for a three-dimensional
    fully nested Archimedean copula of the specified \code{family}.  First, for
    each \eqn{theta_0}{theta0}, \code{n} random variates \eqn{V_0}{V0} are
    generated.  Then, given the particular \eqn{theta_0}{theta0} and the
    realizations \eqn{V_0}{V0}, \code{n} random variates
    \eqn{V_{01}}{V01} are generated for each \eqn{theta_1}{theta1}
    fulfilling the sufficient nesting condition; see \code{paraConstr}
    in \code{\linkS4class{acopula}}.}
  \item{digits}{number of digits for the output.}
  \item{verbose}{logical indicating if \code{nacFrail.time} output
    should generated while the random variates are generated (defaults
    to \code{FALSE}).}
}
\value{
  A \eqn{k \times k}{k x k} matrix of user run time measurements in milliseconds
  (\code{1000*\link{system.time}(.)[1]}) where \eqn{k} is \code{length(taus)}.
  The first column contains the run times for generating the
  \eqn{V_0}{V0}s.  For the submatrix that remains if the first column is
  removed, row \eqn{i} (for \eqn{{theta_0}_i}{theta0[i]}) contains the run times
  for the \eqn{V_{01}}{V01}s for a particular \eqn{theta_0}{theta0} and all the
  admissible \eqn{theta_1}{theta1}s.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
  The class \code{\linkS4class{acopula}} and our predefined \code{"acopula"}
  family objects in \code{\link{acopula-families}}.  For some timings on a
  standard notebook, see \code{\link{demo}(timings)} (or the file
  \file{timings.R} in the demo folder).
}
\examples{
## takes about 7 seconds:% so we rather test a much smaller set in R CMD check
\dontrun{nacFrail.time(10000, "Gumbel", tauso= c(0.05,(1:9)/10, 0.95))}
\dontshow{nacFrail.time(1000, "Gumbel", taus= c(0.5,1,6,9)/10)}
}
\keyword{utilities}
