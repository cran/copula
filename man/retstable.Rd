\name{retstable}
\alias{retstable}
\alias{retstableR}
\title{Sampling Exponentially Tilted Stable Distributions}
\description{
  Generating random variates of an exponentially tilted stable
  distribution of the form
  \deqn{\tilde{S}(\alpha, 1, (\cos(\alpha\pi/2)V_0)^{1/\alpha},
    V_0\mathbf{1}_{\{\alpha=1\}}, h\mathbf{1}_{\{\alpha\ne 1\}}; 1),
  }{tS(alpha, 1, (cos(alpha*pi/2)V0)^(1/alpha), V0*1_(alpha==1),
    h*1_(alpha!=1)),}
  with parameters
  \eqn{\alpha\in(0,1]}{alpha in (0,1]}, \eqn{V_0\in(0,\infty)}{V0 in (0,Inf)},
  and \eqn{h\in[0,\infty)}{h in [0,Inf)} and corresponding
  Laplace-Stieltjes transform
  \deqn{\exp(-V_0((h+t)^\alpha-h^\alpha)),\ t\in[0,\infty];
  }{exp(-V0((h+t)^alpha-h^alpha)), t in [0,Inf];}
  see the references for more details about this distribution.
}
\usage{
retstable(alpha, V0, h = 1, method = NULL)
retstableR(alpha, V0, h = 1)
}
\arguments{
  \item{alpha}{parameter in \eqn{(0,1]}.}
  \item{V0}{vector of values in \eqn{(0,\infty)}{(0,Inf)} (for example, when
    sampling nested Clayton copulas, these are random variates from
    \eqn{F_0}{F0}), that is, the distribution corresponding to
    \eqn{\psi_0}{psi0}.}
  \item{h}{parameter in \eqn{[0,\infty)}{[0,Inf)}.}
  \item{method}{a character string denoting the method to use, currently
    either \code{"MH"} (Marius Hofert's algorithm) or \code{"LD"} (Luc Devroye's
 	 algorithm).  By default, when \code{NULL}, a smart choice is
	 made to use the fastest of these methods depending on the
	 specific values of \eqn{V_0}{V0}.}
}
\details{
  \code{retstableR} is a pure \R version of \code{"MH"}, however, not as fast as
  \code{retstable} (implemented in C, based on both methods) and therefore not
  recommended in simulations when run time matters.
}
\value{
  A vector of variates from \eqn{\tilde{S}(\alpha, 1, .....)}{%
    tS(alpha, 1, .....)}; see above.
}
\author{Marius Hofert, Martin Maechler}
\seealso{
 \code{\link{rstable1}} for sampling stable distributions.
}
\references{
  Devroye, L. (2009)
  Random variate generation for exponentially and polynomially
  tilted stable distributions,
  \emph{ACM Transactions on Modeling and Computer Simulation} \bold{19}, 18, 1--20.

  Hofert, M. (2011)
  Efficiently sampling nested Archimedean copulas,
  \emph{Computational Statistics & Data Analysis} \bold{55}, 57--70.

  Hofert, M. (2012),
  Sampling exponentially tilted stable distributions,
  \emph{ACM Transactions on Modeling and Computer Simulation} \bold{22}, 1.
}
\examples{
## Draw random variates from an exponentially tilted stable distribution
## with given alpha, V0, and h = 1
alpha <- .2
V0 <- rgamma(200, 1)
rETS <- retstable(alpha, V0)

## Distribution plot the random variates -- log-scaled
hist(log(rETS), prob=TRUE)
lines(density(log(rETS)), col=2)
rug (log(rETS))
}
\keyword{distribution}
