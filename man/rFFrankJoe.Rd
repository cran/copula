\name{rFFrankJoe}
\alias{rFFrank}
\alias{rFJoe}
\title{Sampling Distribution F for Frank and Joe}
\description{
  Generate a vector of variates \eqn{V \sim F}{V ~ F} from the distribution
  function \eqn{F} with Laplace-Stieltjes transform
  \deqn{(1-(1-\exp(-t)(1-e^{-\theta_1}))^\alpha)/(1-e^{-\theta_0}),
  }{(1-(1-exp(-t)*(1-e^(-theta1)))^alpha)/(1-e^(-theta0)),}
  for Frank, or
  \deqn{1-(1-\exp(-t))^\alpha,}{1-(1-exp(-t))^alpha} for Joe, respectively,
  where \eqn{\theta_0}{theta0} and \eqn{\theta_1}{theta1} denote two parameters
  of Frank (that is, \eqn{\theta_0,\theta_1\in(0,\infty)}{theta0,theta1
    in (0,Inf)}) and Joe (that is, \eqn{\theta_0,\theta_1\in[1,\infty)}{%
    theta0,theta1 in [1,Inf)}) satisfying
  \eqn{\theta_0\le\theta_1}{theta0 <= theta1}
  and \eqn{\alpha=\theta_0/\theta_1}{alpha=theta0/theta1}.
}
\usage{
rFFrank(n, theta0, theta1, rej)
rFJoe(n, alpha)
}
\arguments{
  \item{n}{number of variates from \eqn{F}.}
  \item{theta0}{parameter \eqn{\theta_0}{theta0}.}
  \item{theta1}{parameter \eqn{\theta_1}{theta1}.}
  \item{rej}{method switch for \code{rFFrank}: if \code{theta0} >
    \code{rej} a rejection from Joe's family (Sibuya distribution) is
    applied (otherwise, a logarithmic envelope is used).}
  \item{alpha}{parameter \eqn{\alpha=
      \theta_0/\theta_1}{alpha = theta0/theta1} in \eqn{(0,1]} for
    \code{rFJoe}.}
}
\value{
  numeric vector of random variates \eqn{V} of length \code{n}.
}
\details{
  \code{rFFrank(n, theta0, theta1, rej)} calls
  \code{\link{rF01Frank}(rep(1,n), theta0, theta1, rej, 1)} and
  \code{rFJoe(n, alpha)} calls \code{\link{rSibuya}(n, alpha)}.
}
\author{Marius Hofert}
\seealso{
  \code{\link{rF01Frank}}, \code{\link{rF01Joe}}, also for references.
  \code{\link{rSibuya}}, and \code{\link{rnacopula}}.
}
\examples{
## Simple definition of the functions:
rFFrank
rFJoe
}
\keyword{distribution}
