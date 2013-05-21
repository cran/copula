\name{Cn}
\alias{Cn}
\alias{C.n}
\alias{dCn}
\title{The Empirical Copula}
\description{
  Given a random sample from a distribution with continuous margins and
  copula C, the empirical copula is a natural nonparametric estimator of
  C.  The function \code{C.n()} computes the empirical copula.

  The function \code{dCn()} approximates first-order partial derivatives
  of the unknown copula.
}
\usage{
 Cn(x, w) ## <-- deprecated!  use  C.n(w, U=pobs(x)) instead!
C.n(u, U, offset=0, method=c("C", "R"))
dCn(u, U, j.ind=1:d, b=1/sqrt(nrow(U)), ...)
}
\arguments{
  \item{u,w}{an \eqn{(m, d)}-\code{\link{matrix}} with elements in
    \eqn{[0,1]} whose lines contain the evaluation points of the
    empirical copula.}
  \item{U,x}{an \eqn{(n, d)}-\code{\link{matrix}} with elements in
    \eqn{[0,1]} and with the same number \eqn{d} of columns as
    \code{u}.  The rows of \code{U} contain the (pseudo-)data on which
    the empirical copula is built.  A multivariate random sample can be
    transformed to an appropriate \code{U} via \code{\link{pobs}()}.}
  \item{j.ind}{\code{\link{integer}} vector of indices \eqn{j} between 1
    and \eqn{d} indicating the dimensions with respect to which
    first-order partial derivatives are approximated.}
  \item{b}{\code{\link{numeric}} giving the bandwidth for approximating
    first-order partial derivatives.}
  \item{offset}{scaling factor (defaults to 0) of the form
    \code{sum(..)/(n+offset)}.}
  \item{method}{\code{\link{character}} string indicating which method
    is applied to compute the empirical copula. \code{method="C"} uses a
    an implmentation in C, \code{method="R"} uses an \R implementation.}
  \item{\dots}{additional argumets passed to \code{C.n()}.}
}
\value{
  \code{C.n()} returns a vector of length \code{m} containing the values
  of the empirical copula of \code{U} at \code{u}.

  \code{dCn()} returns a \eqn{(m, l)}-matrix or an \eqn{m}-vector (for
  \eqn{l=1}; here, \eqn{l} is the length of \code{j.ind}), containing
  the approximated first-order partial derivatives of the unknown
  copula at \code{u}.
}
\note{
  The first version of our empirical copula implementation, \code{Cn()},
  had its two arguments \emph{inverse} compared to \code{C.n()}, and is
  deprecated now. You \bold{must} swap the arguments, and possibly use
  \code{\link{pobs}}, i.e, \bold{ instead of \code{Cn(x, u)}, use
    \code{C.n(u, U=pobs(x))} !}
}
\details{
  There are several asymptotically equivalent definitions of the
  empirical copula.  Here, the empirical copula is simply defined as
  the empirical distribution function computed from the
  pseudo-observations, that is,
  \deqn{C_n(\bm{u})=\frac{1}{n}\sum_{i=1}^n\mathbf{1}_{\{\hat{\bm{U}}_i\le\bm{u}\}},}{(1/n)
  sum(I_{U_i<=u}, i=1, .., n),}
  where \eqn{\hat{\bm{U}}_i}{U_i}, \eqn{i\in\{1,\dots,n\}}{i=1,..,n},
  denote the pseudo-observations (rows in \code{U}) and \eqn{n} the sample size.

  The approximation for the \eqn{j}th partial derivative of the unknown
  copula \eqn{C} is implemented as, for example, in Rémillard and Scaillet
  (2009), % see also Genest et al. (2012)
  and given by
  \deqn{\hat{\dot{C}}_{jn}(\bm{u})=\frac{C_n(u_1,..,u_{j-1},min(u_j+b,1),u_{j+1},..,u_d)-C_n(u_1,..,u_{j-1},max(u_j-b,0),u_{j+1},..,u_d)}{2b},}{%
    hat(C.)[jn](u) = (C[n](u[1], .., u[j-1], min(u[j]+b, 1), u[j+1], .., u[d]) -
                      C[n](u[1], .., u[j-1], max(u[j]-b, 0), u[j+1], .., u[d])) / (2b),}
  where \eqn{b} denotes the bandwidth and \eqn{C_n}{C[n]} the empirical copula.
}
\author{Ivan Kojadinovic, Marius Hofert}
\references{
  R\enc{ü}{u}schendorf, L. (1976).  Asymptotic distributions of
  multivariate rank order statistics,
  \emph{Annals of Statistics} \bold{4}, 912--923.

  Deheuvels, P. (1979).  La fonction de d\enc{é}{e'}pendance empirique et ses
  propri\enc{été}{e'te'}s: un test non param\enc{é}{e'}trique d'ind\enc{é}{e'}pendance,
  \emph{Acad. Roy. Belg. Bull. Cl. Sci.}, 5th Ser. \bold{65}, 274--292.

  Deheuvels, P. (1981).  A non parametric test for independence,
  \emph{Publ. Inst. Statist. Univ. Paris} \bold{26}, 29--50.

  R\enc{é}{e}millard, B. and Scaillet, O. (2009). Testing for equality
  between two copulas. \emph{Journal of Multivariate Analysis}, 100(3),
  pages 377-386.

  %% Genest, C., G. \enc{Nešlehová}{Neslehova}, J., and Quessy, J.-F. (2012).
  %% Tests of symmetry for bivariate copulas,
  %% \emph{Ann. Inst. Stat. Math.} \bold{64}, 811--834.
}
\seealso{\code{\link{pobs}()} for computing pseudo-observations,
  \code{\link{pCopula}()} for evaluating a copula.}
\examples{
n <- 100
d <- 3
family <- "Gumbel"
theta <- 2
cop <- onacopulaL(family, list(theta=theta, 1:d))
set.seed(1)
U <- rCopula(n, cop)

## random points were to evaluate the empirical copula
u <- matrix(runif(n*d), n, d)
ec <- C.n(u, U=U)

## compare with true distribution function
mean(abs(pCopula(u, copula=cop)-ec)) # increase n to decrease this error

## compare the empirical copula and the true copula
## on the diagonal of the unit square
Cn. <- function(x) C.n(do.call(cbind, rep(list(x), d)), U=U)
curve(Cn., 0, 1, main=paste("Diagonal of a", family, "copula"),
      xlab="u", ylab=expression(italic(C)[n](italic(u),..,italic(u))))
pC <- function(x) pCopula(do.call(cbind, rep(list(x), d)), cop)
curve(pC, lty=2, add=TRUE)
legend("topleft", lty=1:2, bty="n", inset=0.02,
       legend=c(expression(italic(C)[n]), expression(italic(C))))

## check the empirical copula with its Kendall distribution function
plot( pK(C.n(U, U=U), cop=cop@copula, d=d) ) # must be uniform

## approximate partial derivatives w.r.t. the 2nd and 3rd component
j.ind <- 2:3
der23 <- dCn(u, U=pobs(U), j.ind=j.ind)
der23. <- copula:::dCdu(archmCopula(family, param=theta, dim=d), u=u)[,j.ind]
summary(as.vector(abs(der23-der23.))) # approximation error summary
\dontshow{% check backward compatibility
  U <- U[1:64 ,]# such that m != n
  stopifnot(suppressWarnings( ## deprecation warning ..
    identical(C.n(u, pobs(U)),
              Cn (U, u))))
}
}
\keyword{multivariate}

