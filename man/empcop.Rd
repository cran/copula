\name{C.n}
\title{The Empirical Copula}
\alias{C.n}
\alias{Cn}% <- deprecated, but of course keep here
\alias{F.n}
%% \alias{.Fn}
\alias{dCn}
\description{
  Given pseudo-observations from a distribution with continuous margins and
  copula C, the empirical copula is the empirical distribution function
  of these pseudo-observations.  It is thus a natural nonparametric estimator of
  C.  The function \code{C.n()} computes the empirical copula.

  The function \code{dCn()} approximates first-order partial derivatives
  of the unknown copula.

  The function \code{F.n()} computes the empirical distribution function
  of a multivariate sample.  Note that \code{C.n(u, X, *)} simply calls
  \code{F.n(u, pobs(X), *)} after checking \code{u}.
}
\usage{
C.n(u, X, offset=0, method=c("C", "R"))
dCn(u, U, j.ind=1:d, b=1/sqrt(nrow(U)), ...)

F.n(x, X, offset=0, method=c("C", "R"))

Cn(x, w) ## <-- deprecated!  use  C.n(w, x) instead!
}
\arguments{
  \item{u,w}{an \eqn{(m, d)}-\code{\link{matrix}} with elements in
    \eqn{[0,1]} whose rows contain the evaluation points of the
    empirical copula.}
  \item{x}{an \eqn{(m, d)}-\code{\link{matrix}} whose rows
    contain the evaluation points of the empirical distribution function.}
  \item{U}{for \code{dCN()} only: an \eqn{(n,d)}-\code{\link{matrix}} with
    elements in \eqn{[0,1]} and with the same number \eqn{d} of columns as
    \code{u}.  The rows of \code{U} are the pseudo-observations based
    on which the empirical copula is computed.}
  \item{X}{(and \code{x} and \code{U} for \code{Cn()}:) an \eqn{(n,
      d)}-\code{\link{matrix}} with the same number \eqn{d} of columns as
    \code{x}.  Recall that a multivariate random sample \code{X} can be
    transformed to an appropriate \code{U} via \code{\link{pobs}()}.}
  \item{j.ind}{\code{\link{integer}} vector of indices \eqn{j} between 1
    and \eqn{d} indicating the dimensions with respect to which
    first-order partial derivatives are approximated.}
  \item{b}{\code{\link{numeric}} giving the bandwidth for approximating
    first-order partial derivatives.}
  \item{offset}{used in scaling the result which is of the form
    \code{sum(....)/(n+offset)}; defaults to zero.}
  \item{method}{\code{\link{character}} string indicating which method
    is applied to compute the empirical cumulative distribution function
    or the empirical copula.
    \code{method="C"} uses a an implementation in C,
    \code{method="R"} uses a pure \R implementation.}
  \item{\dots}{additional arguments passed to \code{dCn()}.}
}
\value{
  \code{C.n()} returns the empirical copula at \code{u} (that is, the
  empirical distribution function of the observations \code{U}
  evaluated at \code{u}).  The \emph{name} \dQuote{U} is a slight
  misnomer here, for back compatibility reasons, as \code{U} typically
  contains the original observations \code{X} here.  \code{F.n()}
  returns the empirical distribution function of \code{X} evaluated at
  \code{x}.
  %% The result is slightly wrapped as \code{"mvFn"} (S3) class object,
  %% with a \code{\link{persp}} and \code{\link{print}} method, but behaves
  %% otherwise a simple \code{\link{numeric}} vector.
  %% As pure computational tool, \code{.Fn()} returns the same raw numeric
  %% vector.

  \code{dCn()} returns a vector (\code{length(j.ind)} is 1) or a matrix
  (with number of columns equal to \code{length(j.ind)}), containing
  the approximated first-order partial derivatives of the unknown
  copula at \code{u} with respect to the arguments in \code{j.ind}.
}
\note{
  The first version of our empirical copula implementation, \code{Cn()},
  had its two arguments \emph{reversed} compared to \code{C.n()}, and is
  deprecated now.  You \bold{must} swap its arguments to transform to
  new code.
}
\details{
  There are several asymptotically equivalent definitions of the
  empirical copula. As mentioned above, the empirical copula
  \code{C.n()} is simply defined as the empirical distribution function
  computed from the pseudo-observations, that is,
  \deqn{C_n(\bm{u})=\frac{1}{n}\sum_{i=1}^n\mathbf{1}_{\{\hat{\bm{U}}_i\le\bm{u}\}},}{(1/n)
  sum(I_{U_i<=u}, i=1, .., n),}
  where \eqn{\hat{\bm{U}}_i}{U_i}, \eqn{i\in\{1,\dots,n\}}{i=1,..,n},
  denote the pseudo-observations (rows in \code{U}) and \eqn{n} the
  sample size. Internally, \code{C.n()} is just a wrapper for
  \code{F.n()} and is expected to be fed with the pseudo-observations.

  The approximation for the \eqn{j}th partial derivative of the unknown
  copula \eqn{C} is implemented as, for example, in R\enc{é}{e'}millard and Scaillet
  (2009), % see also Genest et al. (2012)
  and given by
  \deqn{\hat{\dot{C}}_{jn}(\bm{u})=\frac{C_n(u_1,..,u_{j-1},min(u_j+b,1),u_{j+1},..,u_d)-C_n(u_1,..,u_{j-1},max(u_j-b,0),u_{j+1},..,u_d)}{2b},}{%
    hat(C.)[jn](u) = (C[n](u[1], .., u[j-1], min(u[j]+b, 1), u[j+1], .., u[d]) -
                      C[n](u[1], .., u[j-1], max(u[j]-b, 0), u[j+1], .., u[d])) / (2b),}
  where \eqn{b} denotes the bandwidth and \eqn{C_n}{C[n]} the empirical copula.
}
%\author{Ivan Kojadinovic and Marius Hofert}
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
  \code{\link{pCopula}()} for evaluating a copula.
}
\examples{
## Generate data X (from a meta-Gumbel model with N(0,1) margins)
n <- 100
d <- 3
family <- "Gumbel"
theta <- 2
cop <- onacopulaL(family, list(theta=theta, 1:d))
set.seed(1)
X <- qnorm(rCopula(n, cop)) # meta-Gumbel data with N(0,1) margins

## Random points were to evaluate the empirical copula
u <- matrix(runif(n*d), n, d)
ec <- C.n(u, X)

## Compare the empirical copula with the true copula
mean(abs(pCopula(u, copula=cop)-ec)) # ~= 0.012 -- increase n to decrease this error

## Compare the empirical copula with F.n(pobs())
U <- pobs(X) # pseudo-observations
stopifnot(identical(ec, F.n(u, X=pobs(U)))) # even identical

## Compare the empirical copula based on U at U with the Kendall distribution
## Note: Theoretically, C(U) ~ K, so K(C_n(U, U=U)) should approximately be U(0,1)
plot(pK(C.n(U, X), cop=cop@copula, d=d))

## Compare the empirical copula and the true copula on the diagonal
C.n.diag <- function(u) C.n(do.call(cbind, rep(list(u), d)), X=X) # diagonal of C_n
C.diag <- function(u) pCopula(do.call(cbind, rep(list(u), d)), cop) # diagonal of C
curve(C.n.diag, from=0, to=1, # empirical copula diagonal
      main=paste("True vs empirical diagonal of a", family, "copula"),
      xlab="u", ylab=expression("True C(u,..,u) and empirical"~C[n](u,..,u)))
curve(C.diag, lty=2, add=TRUE) # add true copula diagonal
legend("bottomright", lty=2:1, bty="n", inset=0.02,
       legend=c("C", expression(C[n])))

## Approximate partial derivatives w.r.t. the 2nd and 3rd component
j.ind <- 2:3 # indices w.r.t. which the partial derivatives are computed
## Partial derivatives based on the empirical copula and the true copula
der23 <- dCn(u, U=pobs(U), j.ind=j.ind)
der23. <- copula:::dCdu(archmCopula(family, param=theta, dim=d), u=u)[,j.ind]
## Approximation error
summary(as.vector(abs(der23-der23.)))
\dontshow{% check backward compatibility
  U <- U[1:64 ,]# such that m != n
  stopifnot(suppressWarnings( ## deprecation warning ..
    identical(C.n(u, U),
              Cn (U, u))))
}
## For an example of using F.n(), see help(mvdc)% ./Mvdc.Rd
}
\keyword{multivariate}
