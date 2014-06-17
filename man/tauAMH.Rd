\name{tauAMH}
\alias{tauAMH}
\alias{tauJoe}
\title{Ali-Mikhail-Haq ("AMH")'s and Joe's Kendall's Tau}
\description{
  Compute Kendall's Tau of an Ali-Mikhail-Haq ("AMH") or Joe Archimedean
  copula with parameter \code{theta}.  In both cases, analytical
  expressions are available, but need alternatives in some cases.
  \describe{
    \item{\code{tauAMH()}:}{Analytically, given as
      \deqn{1-\frac{2((1-\theta)^2\log(1-\theta) + \theta)}{3\theta^2},}{
	1 - 2((1-t)(1-t)log(1-t) + t)/(3*t^2),}
      for \code{theta}\eqn{=\theta}{=t};
      numerically, care has to be taken when \eqn{\theta \to 0}{t -> 0},
      avoiding accuracy loss already, for example, for \eqn{\theta}{t} as
      large as \code{theta = 0.001}.}
    \item{\code{tauJoe()}:}{
      Analytically,
      \deqn{1- 4\sum_{k=1}^\infty\frac{1}{k(\theta k+2)(\theta(k-1)+2)},}{%
	1- 4 sum{k=1:Inf; 1/(k(t*k+2)(t(k-1)+2))},}
      the infinite sum can be expressed by three \eqn{\psi()}
      (\code{\link{psigamma}}) function terms.
    }
  }
}
\usage{
tauAMH(theta)
tauJoe(theta, method = c("hybrid", "digamma", "sum"), noTerms=446)
}
\arguments{
  \item{theta}{numeric vector with values in \eqn{[-1,1]} for AMH, or
    \eqn{[0.238734, Inf)} for Joe.}
  \item{method}{string specifying the method for \code{tauJoe()}.  Use the
    default, unless for research about the method.  Up to \pkg{copula}
    version 0.999-0, the only (implicit) method was \code{"sum"}.}
  \item{noTerms}{the number of summation terms for the \code{"sum"}
    method; its default, \code{446} gives an absolute error smaller
    than \eqn{10^{-5}}.}
}
\value{
  a vector of the same length as \code{theta} (\eqn{= \theta}), with
  \eqn{\tau} values

  for \code{tauAMH}: in \eqn{[(5 - 8 log 2)/3, 1/3] ~= [-0.1817, 0.3333]},
  of \eqn{\tau_A(\theta) = 1 -
    2(\theta+(1-\theta)^2\log(1-\theta))/(3\theta^2)}{%
    tau.A(t) = 1 - 2*((1-t)*(1-t)*log(1-t) + t) / (3*t^2)},
  numerically accurately, to at least around 12 decimal digits.

  for  \code{tauJoe}: in [-1,1].
}
\details{
\describe{
  \item{\code{tauAMH()}:}{
    For small \code{theta} (\eqn{=\theta}), we use Taylor series
    approximations of up to order 7,
    \deqn{\tau_A(\theta) = \frac{2}{9}\theta\Bigl(1 + \theta\Bigl(\frac 1 4 +
      \frac{\theta}{10}\Bigl(1 + \theta\Bigl(\frac 1 2 + \theta \frac 2
      7\Bigr) \Bigr)\Bigr)\Bigr) + O(\theta^6),}{%
      tau[A](th) = 2*th/9 *(1 + th*(1/4 + th/10* (1 + th*(1/2 + th*2/7))))
      + O(th^6),}
    where we found that dropping the last two terms (e.g., only using 5 terms
    from the \eqn{k=7} term Taylor polynomial) is actually numerically
    advantageous.}
  \item{\code{tauJoe()}:}{
    The \code{"sum"} method simply replaces the infinite sum by a finite
    sum (with \code{noTerms} terms.  The more accurate or faster methods,
    use analytical summation formulas, using the \code{\link{digamma}}
    aka \eqn{\psi} function, see, e.g.,
    \url{http://en.wikipedia.org/wiki/Digamma_function#Series_formula}.
%% FIXME: Martin, use your notes (Aug.7, 2012; oranges Heft) to write a
%%      small "paper" / vignette, where everything is nicely derived,
%% notably the cases theta=2, and theta ~= 2.

    The smallest sensible \eqn{\theta} value, i.e., \code{th} for which
    \code{tauJoe(th) == -1} is easily determined via
    \code{str(uniroot(function(th) tauJoe(th)-(-1), c(0.1, 0.3), tol = 1e-17), digits=12)}
    to be \code{0.2387339899}.}
}}
\author{Martin Maechler}% copAMH: June 2010;  copJoe: August 2012
\seealso{
  \code{\link{acopula-families}}, and their class definition,
  \code{"\linkS4class{acopula}"}.
}
\examples{
tauAMH(c(0, 2^-40, 2^-20))
curve(tauAMH,  0, 1)
curve(tauAMH, -1, 1)# negative taus as well
curve(tauAMH, 1e-12, 1, log="xy") # linear, tau ~= 2/9*theta in the limit

curve(tauJoe, 1,      10)
curve(tauJoe, 0.2387, 10)# negative taus (*not* valid for Joe: no 2-monotone psi()!)
}
\keyword{distribution}
