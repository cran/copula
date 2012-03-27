\name{math-fun}
\alias{A..Z}
\alias{sinc}
\title{Sinc, Zolotarev's, and Other Mathematical Utility Functions}
\description{
  \code{sinc(x)} computes the \emph{sinc function}
  \eqn{s(x)=\sin(x)/x}{s(x) = sin(x)/x} for \eqn{x\ne 0}{x != 0} and \eqn{s(0) = 1},
  such that \eqn{s()} is continuous, also at \eqn{x = 0}.

  \code{A..Z(x, a)} computes Zolotarev's function to
  the power \code{1-a}.

}
\usage{
sinc(x)
A..Z(x, alpha, I.alpha = 1 - alpha)
}
\arguments{
  \item{x}{\code{\link{numeric}} argument in \eqn{[0,\pi]}{[0,pi]}, typically a vector.}
  \item{alpha}{parameter in (0,1].}
  \item{I.alpha}{must be \code{ = 1 - alpha}, maybe more accurately
    when \code{alpha} is very close to 1.}
}
\details{
  For more details about Zolotarev's function, see, for example, Devroye (2009).
}
\value{
  \code{A..Z(x,alpha)} is \eqn{\tilde A_{Z}(x,\alpha)}{A~Z(x,alpha)},
  defined as
  \deqn{\frac{\sin(\alpha x)^\alpha\sin((1-\alpha)x)^{1-\alpha}}{\sin(x)},\
    x\in[0,\pi],}{%
    sin(alpha*x)^alpha * sin((1-alpha)*x)^(1-alpha) / sin(x),  x in [0,pi],}

  where \eqn{\alpha\in(0,1]}{alpha in (0,1]} is \code{alpha}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{retstable}} internally makes use of these functions.
}
\references{
  Devroye, L. (2009)
  Random variate generation for exponentially and polynomially tilted
  stable distributions,
  \emph{ACM Transactions on Modeling and Computer Simulation} \bold{19}, 18, 1--20.
}
\examples{
curve(sinc, -15,25); abline(h=0,v=0, lty=2)
curve(A..Z(x, 0.25), xlim = c(-4,4),
      main = "Zolotarev's function A(x) ^ 1-alpha")
}
\keyword{math}
