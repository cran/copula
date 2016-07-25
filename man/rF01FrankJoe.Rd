\name{rF01FrankJoe}
\alias{rF01Frank}
\alias{rF01Joe}
\title{Sample Univariate Distributions Involved in Nested Frank and Joe Copulas}
\description{
  \code{rF01Frank}:
  Generate a vector of random variates \eqn{V_{01}\sim F_{01}}{V01 ~ F01}
  with Laplace-Stieltjes transform
  \deqn{\psi_{01}(t;V_0)=
    \Bigl(\frac{1-(1-\exp(-t)(1-e^{-\theta_1}))^{\theta_0/\theta_1}}{%
      1-e^{-\theta_0}}\Bigr)^{V_0}.}{%
    psi01(t;V0) = ((1-(1-exp(-t)*(1-e^-theta1))^(theta0/theta1))/
    (1-e^-theta0))^V0.
  }
  for the given	realizations \eqn{V_0}{V0} of Frank's \eqn{F_0}{F0} and the
  parameters \eqn{\theta_0,\theta_1\in(0,\infty)}{theta0, theta1 in (0,Inf)}
  such that \eqn{\theta_0\le\theta_1}{theta0 <= theta1}.  This distribution
  appears on sampling nested Frank copulas.  The parameter \code{rej} is used
  to determine the cut-off point of two algorithms that are involved in
  sampling \eqn{F_{01}}{F01}.  If
  \eqn{\code{rej} < V_0\theta_0(1-e^{-\theta_0})^{V_0-1}}{%
    rej < V0*theta_0*(1-e^{-theta0})^(V0-1)}
  a rejection from \eqn{F_{01}}{F01} of Joe is applied (see
  \code{rF01Joe}; the meaning of the parameter \code{approx} is
  explained below), otherwise a sum is sampled with a logarithmic
  envelope for each summand.

  \code{rF01Joe}:
  Generate a vector of random variates \eqn{V_{01}\sim F_{01}}{V01 ~ F01}
  with Laplace-Stieltjes transform
  \deqn{\psi_{01}(t;V_0)=(1-(1-\exp(-t))^\alpha)^{V_0}.}{%
    psi01(t;V0) = (1-(1-exp(-t))^alpha)^V0.}
  for the given realizations \eqn{V_0}{V0} of Joe's \eqn{F_0}{F0} and
  the parameter \eqn{\alpha\in(0,1]}{alpha in (0,1]}.  This distribution
  appears on sampling nested Joe copulas.  Here,
  \eqn{\alpha=\theta_0/\theta_1}{alpha = theta0/theta1}, where
  \eqn{\theta_0,\theta_1\in[1,\infty)}{theta0, theta1 in [1,Inf)} such
  that \eqn{\theta_0\le\theta_1}{theta0 <= theta1}.  The parameter
  \code{approx} denotes the largest number of summands in the
  sum-representation of \eqn{V_{01}}{V01} before the asymptotic
  \deqn{V_{01}=V_0^{1/\alpha}S(\alpha,1,\cos^{1/\alpha}(\alpha\pi/2),
    \mathbf{1}_{\{\alpha=1\}};1)}{%
    V01 = V0^(1/alpha) S(alpha,1, cos^(1/alpha)(alpha*pi/2), 1_(alpha==1); 1)}
  is used to sample \eqn{V_{01}}{V01}.
}
\usage{
rF01Frank(V0, theta0, theta1, rej, approx)
rF01Joe(V0, alpha, approx)
}
\arguments{
  \item{V0}{a vector of random variates from \eqn{F_0}{F0}.}
  \item{theta0, theta1, alpha}{parameters
    \eqn{\theta_0,\theta_1}{theta0, theta1}
    and \eqn{\alpha}{alpha} as described above.}
  \item{rej}{parameter value as described above.}
  \item{approx}{parameter value as described above.}
}
\value{
  A vector of positive \code{\link{integer}}s of length \code{n}
  containing the generated random variates.
}
%\author{Marius Hofert and Martin Maechler}
\references{
  Hofert, M. (2011).
  Efficiently sampling nested Archimedean copulas.
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.
}
\seealso{
  \code{\link{rFFrank}}, \code{\link{rFJoe}}, \code{\link{rSibuya}},
  and \code{\link{rnacopula}}.

  \code{\link{rnacopula}}
}
\examples{
## Sample n random variates V0 ~ F0 for Frank and Joe with parameter
## chosen such that Kendall's tau equals 0.2 and plot histogram
n <- 1000
theta0.F <- copFrank@iTau(0.2)
V0.F <- copFrank@V0(n,theta0.F)
hist(log(V0.F), prob=TRUE); lines(density(log(V0.F)), col=2, lwd=2)
theta0.J <- copJoe@iTau(0.2)
V0.J <- copJoe@V0(n,theta0.J)
hist(log(V0.J), prob=TRUE); lines(density(log(V0.J)), col=2, lwd=2)

## Sample corresponding V01 ~ F01 for Frank and Joe and plot histogram
## copFrank@V01 calls rF01Frank(V0, theta0, theta1, rej=1, approx=10000)
## copJoe@V01 calls rF01Joe(V0, alpha, approx=10000)
theta1.F <- copFrank@iTau(0.5)
V01.F <- copFrank@V01(V0.F,theta0.F,theta1.F)
hist(log(V01.F), prob=TRUE); lines(density(log(V01.F)), col=2, lwd=2)
theta1.J <- copJoe@iTau(0.5)
V01.J <- copJoe@V01(V0.J,theta0.J,theta1.J)
hist(log(V01.J), prob=TRUE); lines(density(log(V01.J)), col=2, lwd=2)
}
\keyword{distribution}
