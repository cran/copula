\name{empCopula}
\title{The Empirical Copula}
\alias{F.n}
\alias{dCn}
\alias{C.n}
\alias{empCopula}
\alias{Cn}% <- deprecated, but of course kept here
\alias{toEmpMargins}
\description{
  Computes the empirical copula (according to a provided method) and
  auxiliary tools.
}
\usage{
empCopula(X, smoothing = c("none", "beta", "checkerboard"), offset = 0,
          ties.method = c("max", "average", "first", "last", "random", "min"))
C.n(u, X, smoothing = c("none", "beta", "checkerboard"), offset = 0,
    ties.method = c("max", "average", "first", "last", "random", "min"))
dCn(u, U, j.ind = 1:d, b = 1/sqrt(nrow(U)), ...)
F.n(x, X, offset = 0, smoothing = c("none", "beta", "checkerboard"))
Cn(x, w) ## <-- deprecated!  use  C.n(w, x) instead!
toEmpMargins(U, x, ...)
}
\arguments{
  \item{X}{an \eqn{(n, d)}-\code{\link{matrix}} of pseudo-observations
    with \eqn{d} columns (as \code{x} or \code{u}). Recall that a
    multivariate random sample can be transformed to pseudo-observations
    via \code{\link{pobs}()}. For \code{F.n()} and if \code{smoothing !=
      "none"}, \code{X} can also be a general, multivariate sample,
    in which case the empirical distribution function is computed.}
  \item{u,w}{an \eqn{(m, d)}-\code{\link{matrix}} with elements in
    \eqn{[0,1]} whose rows contain the evaluation points of the
    empirical copula.}
  \item{U}{an \eqn{(n,d)}-\code{\link{matrix}} of
    pseudo- (or copula-)observations (elements in \eqn{[0,1]},
    same number \eqn{d} of columns as \code{u} (for \code{dCn()}))
    or \code{x} (for \code{toEmpMargins()}).}
  \item{x}{an \eqn{(m, d)}-\code{\link{matrix}} whose rows
    contain the evaluation points of the empirical distribution< function
    (if \code{smoothing = "none"}) or copula (if \code{smoothing != "none"}).}
  \item{smoothing}{\code{\link{character}} string specifying whether
    the empirical distribution function (for \code{F.n()}) or copula
    (for \code{C.n()}) is computed (if \code{smoothing = "none"}),
    or whether the empirical beta copula (\code{smoothing = "beta"})
    or the empirical checkerboard  copula (\code{smoothing =
      "checkerboard"}) is computed.}
  \item{ties.method}{\code{\link{character}} string specifying how ranks should be computed if
    there are ties in any of the coordinate samples of \code{x}; passed
    to \code{\link{pobs}}.}
  \item{j.ind}{\code{\link{integer}} vector of indices \eqn{j} between 1
    and \eqn{d} indicating the dimensions with respect to which
    first-order partial derivatives are approximated.}
  \item{b}{\code{\link{numeric}} giving the bandwidth for approximating
    first-order partial derivatives.}
  \item{offset}{used in scaling the result which is of the form
    \code{sum(....)/(n+offset)}; defaults to zero.}
  \item{\dots}{additional arguments passed to \code{dCn()} or
    \code{sort()} underlying \code{toEmpMargins()}.}
}
\value{
  \code{empCopula()} is the constructor for objects of class
  \code{\linkS4class{empCopula}}.

  \code{C.n()} returns the empirical copula of the pseudo-observations
  \code{X} evaluated at \code{u} (or a smoothed version of it).
  %(that is, the
  %empirical distribution function of the observations \code{U}
  %evaluated at \code{u})
  %The \emph{name} \dQuote{U} is a slight
  %misnomer here, for back compatibility reasons, as \code{U} typically
  %contains the original observations \code{X} here.

  \code{dCn()} returns a vector (\code{length(j.ind)} is 1) or a matrix
  (with number of columns equal to \code{length(j.ind)}), containing
  the approximated first-order partial derivatives of the unknown
  copula at \code{u} with respect to the arguments in \code{j.ind}.

  \code{F.n()} returns the empirical distribution function of \code{X}
  evaluated at \code{x} if \code{smoothing = "none"}, the empirical beta
  copula evaluated at \code{x} if \code{smoothing = "beta"} and the empirical
  checkerboard copula evaluated at \code{x} if \code{smoothing =
    "checkerboard"}.
  %% The result is slightly wrapped as \code{"mvFn"} (S3) class object,
  %% with a \code{\link{persp}} and \code{\link{print}} method, but behaves
  %% otherwise a simple \code{\link{numeric}} vector.
  %% As pure computational tool, \code{.Fn()} returns the same raw numeric
  %% vector.

  \code{toEmpMargins()} transforms the copula sample \code{U} to the
  empirical margins based on the sample \code{x}.
}
\note{
  The first version of our empirical copula implementation, \code{Cn()},
  had its two arguments \emph{reversed} compared to \code{C.n()}, and is
  deprecated now.  You \bold{must} swap its arguments to transform to
  new code.

  The use of the two smoothed versions assumes implicitly no
  ties in the component samples of the data.
}
\details{
  Given pseudo-observations from a distribution with continuous margins
  and copula C, the \emph{empirical copula} is the (default) empirical
  distribution function of these pseudo-observations.  It is thus a
  natural nonparametric estimator of C. The function \code{C.n()}
  computes the empirical copula or two alternative smoothed versions of
  it: the \emph{empirical beta copula} or the \emph{empirical
  checkerboard copula}; see Eqs. (2.1) and (4.1) in Segers, Sibuya and
  Tsukahara (2017), and the references therein. \code{empCopula()}
  is the constructor of an object of class \code{\linkS4class{empCopula}}.

  The function \code{dCn()} approximates first-order partial derivatives
  of the unknown copula using the empirical copula.

  The function \code{F.n()} computes the empirical distribution function
  of a multivariate sample. Note that \code{C.n(u, X, smoothing="none",
  *)} simply calls \code{F.n(u, pobs(X), *)} after checking \code{u}.

  There are several asymptotically equivalent definitions of the
  empirical copula. \code{C.n(, smoothing = "none")} is simply defined
  as the empirical distribution function computed from the
  pseudo-observations, that is,
  \deqn{C_n(\bm{u})=\frac{1}{n}\sum_{i=1}^n\mathbf{1}_{\{\hat{\bm{U}}_i\le\bm{u}\}},}{(1/n)
  sum(I_{U_i<=u}, i=1, .., n),} where \eqn{\hat{\bm{U}}_i}{U_i},
  \eqn{i\in\{1,\dots,n\}}{i=1,..,n}, denote the pseudo-observations and
  \eqn{n} the sample size. Internally, \code{C.n(,smoothing = "none")}
  is just a wrapper for \code{F.n()} and is expected to be fed with the
  pseudo-observations.

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

  Segers, J., Sibuya, M. and Tsukahara, H. (2017). The Empirical Beta
  Copula. \emph{Journal of Multivariate Analysis}, 155, pages 35--51,
  \url{http://arxiv.org/abs/1607.04430}.

  %% Genest, C., G. \enc{Nešlehová}{Neslehova}, J., and Quessy, J.-F. (2012).
  %% Tests of symmetry for bivariate copulas,
  %% \emph{Ann. Inst. Stat. Math.} \bold{64}, 811--834.
}
\seealso{
  \code{\link{pobs}()} for computing pseudo-observations.
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

## Evaluate empirical copula
u <- matrix(runif(n*d), n, d) # random points were to evaluate the empirical copula
ec <- C.n(u, X = X)

## Compare the empirical copula with the true copula
pc <- pCopula(u, copula = cop)
mean(abs(pc - ec)) # ~= 0.012 -- increase n to decrease this error

## The same for the two smoothed versions
beta <- C.n(u, X, smoothing = "beta")
mean(abs(pc - beta))
check <- C.n(u, X, smoothing = "checkerboard")
mean(abs(pc - check))

## Compare the empirical copula with F.n(pobs())
U <- pobs(X) # pseudo-observations
stopifnot(identical(ec, F.n(u, X = pobs(U)))) # even identical

## Compare the empirical copula based on U at U with the Kendall distribution
## Note: Theoretically, C(U) ~ K, so K(C_n(U, U = U)) should approximately be U(0,1)
plot(ecdf(pK(C.n(U, X), cop = cop@copula, d = d)), asp = 1, xaxs="i", yaxs="i")
segments(0,0, 1,1, col=adjustcolor("blue",1/3), lwd=5, lty = 2)
abline(v=0:1, col="gray70", lty = 2)

## Compare the empirical copula and the true copula on the diagonal
C.n.diag <- function(u) C.n(do.call(cbind, rep(list(u), d)), X = X) # diagonal of C_n
C.diag <- function(u) pCopula(do.call(cbind, rep(list(u), d)), cop) # diagonal of C
curve(C.n.diag, from = 0, to = 1, # empirical copula diagonal
      main = paste("True vs empirical diagonal of a", family, "copula"),
      xlab = "u", ylab = quote("True C(u,..,u) and empirical"~C[n](u,..,u)))
curve(C.diag, lty = 2, add = TRUE) # add true copula diagonal
legend("bottomright", lty = 2:1, bty = "n", inset = 0.02,
       legend = expression(C, C[n]))

## Approximate partial derivatives w.r.t. the 2nd and 3rd component
j.ind <- 2:3 # indices w.r.t. which the partial derivatives are computed
## Partial derivatives based on the empirical copula and the true copula
der23 <- dCn(u, U = pobs(U), j.ind = j.ind)
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

## Generate a bivariate empirical copula object (various smoothing methods)
n <- 10 # sample size
d <- 2 # dimension
set.seed(271)
X <- rCopula(n, copula = claytonCopula(3, dim = d))
ecop.orig  <- empCopula(X) # smoothing = "none"
ecop.beta  <- empCopula(X, smoothing = "beta")
ecop.check <- empCopula(X, smoothing = "checkerboard")

## Plot
wireframe2(ecop.orig,  FUN = pCopula, draw.4.pCoplines = FALSE)
wireframe2(ecop.beta,  FUN = pCopula)
wireframe2(ecop.check, FUN = pCopula)
## Density (only exists when smoothing = "beta")
wireframe2(ecop.beta,  FUN = dCopula)

## Transform a copula sample to empirical margins
set.seed(271)
X <- qexp(rCopula(1000, copula = claytonCopula(2))) # multivariate distribution
U <- rCopula(917, copula = gumbelCopula(2)) # new copula sample
X. <- toEmpMargins(U, x = X) # tranform U to the empirical margins of X
plot(X.) # Gumbel sample with empirical margins of X
}
\keyword{multivariate}
