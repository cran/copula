\name{tauAMH}
\alias{tauAMH}
\title{Ali-Mikhail-Haq ("AMH")'s Kendall's tau}
\description{
  Compute Kendall's Tau of an Ali-Mikhail-Haq ("AMH") Archimedean copula
  with parameter \code{theta}.  While that's analytically given
  explicitly, as
  \deqn{1-\frac{2((1-\theta)^2\log(1-\theta) + \theta)}{3\theta^2},}{
    1 - 2((1-t)(1-t)log(1-t) + t)/(3*t^2),}
  for \code{th}\eqn{=\theta}{=t};
  numerically, care has to be taken when \eqn{\theta \to 0}, avoiding accuracy
  loss already, for example, for \eqn{\theta} as large as \code{theta = 0.001}.
}
\details{
  For small \code{th} (\eqn{=\theta}), we use Taylor series
  approximations of up to order 7,
  \deqn{\tau_A(\theta) = \frac{2}{9}\theta\Bigl(1 + \theta\Bigl(\frac 1 4 +
    \frac{\theta}{10}\Bigl(1 + \theta\Bigl(\frac 1 2 + \theta \frac 2
    7\Bigr) \Bigr)\Bigr)\Bigr) + O(\theta^6),}{%
    tau[A](th) = 2*th/9 *(1 + th*(1/4 + th/10* (1 + th*(1/2 + th*2/7))))
  + O(th^6),}
  where found that dropping the last two terms (e.g., only using 5 terms
  from the \eqn{k=7} term Taylor polynomial) is actually numerically
  advantageous.
}
\usage{
tauAMH(th)
}
\arguments{
  \item{th}{numeric vector with values in \eqn{[0,1]}.}
}
\value{
  a vector of the same length as \code{th} (\eqn{= \theta}), with
  \eqn{\tau} values in \eqn{[0, 1/3]}, of
  \eqn{\tau_A(\theta) = 1 - 2(\theta+(1-\theta)^2\log(1-\theta))/(3\theta^2)}{%
    tau.A(th) = 1 - 2*((1-th)*(1-th)*log(1-th) + th) / (3*th^2)}, numerically
  accurately, to at least around 12 decimal digits.
}
\author{Martin Maechler, June 2010}
\seealso{
  \code{\link{acopula-families}}, and their class definition,
  \code{"\linkS4class{acopula}"}.
}
\examples{
tauAMH(c(0, 2^-40, 2^-20))
curve(tauAMH, 0, 1)
curve(tauAMH, 1e-12, 1, log="xy") # linear, tau ~= 2/9*theta in the limit
}
\keyword{distribution}
