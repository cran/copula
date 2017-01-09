\name{K}
\title{Kendall Distribution Function for Archimedean Copulas}
\alias{Kn}
\alias{pK}
\alias{qK}
\alias{dK}
\alias{rK}
\alias{K}% Defunct
\description{
  The \emph{Kendall distribution} of an Archimedean
  copula is defined by
  \deqn{K(u) = P(C(U_1,U_2,\dots,U_d) \le u),}{K(u) = P(C(U[1], U[2], \dots, U[d]) <= u),}
  where \eqn{u \in [0,1]}{u in [0,1]}, and the \eqn{d}-dimensional
  \eqn{(U_1,U_2,\dots,U_d)}{(U[1],U[2],\dots,U[d])} is distributed according
  to the copula \eqn{C}.  Note that the random variable
  \eqn{C(U_1,U_2,\dots,U_d)}{C(U[1],U[2], \dots, U[d])} is known as
  \dQuote{probability integral transform}.  Its distribution function
  \eqn{K} is equal to the identity if \eqn{d = 1}, but is non-trivial for
  \eqn{d \ge 2}{d >= 2}.

  \code{Kn()} computes the empirical Kendall distribution function,
  \code{pK()} the distribution function (so \eqn{K()} itself),
  \code{qK()} the quantile function, \code{dK()} the density, and
  \code{rK()} random number generation from \eqn{K()} for an Archimedean
  copula.
}
\usage{
% K(u, copula, d, n.MC=0, log.p=FALSE) # is deprecated
Kn(u, x) # empirical Kendall distribution function
dK(u, copula, d, n.MC = 0, log.p = FALSE) # density
pK(u, copula, d, n.MC = 0, log.p = FALSE) # df
qK(p, copula, d, n.MC = 0, log.p = FALSE, # quantile function
   method = c("default", "simple", "sort", "discrete", "monoH.FC"),
   u.grid, \dots)
rK(n, copula, d) # random number generation
}
\arguments{
  \item{u}{evaluation point(s) (in \eqn{[0,1]}).}
  \item{x}{data (in the \eqn{d}-dimensional space) based on which the
    Kendall distribution function is estimated.}
  \item{copula}{\code{\linkS4class{acopula}} with specified parameter, or
    (currently for \code{rK} only) a \code{\linkS4class{outer_nacopula}}.}
  \item{d}{dimension (not used when \code{copula} is an
    \code{\linkS4class{outer_nacopula}}).}
  \item{n.MC}{\code{\link{integer}}, if positive, a Monte Carlo approach
    is applied with sample size equal to \code{n.MC} to evaluate the
    generator derivatives involved; otherwise
    (\code{n.MC = 0}) the exact formula is used based on the generator
    derivatives as found by Hofert et al. (2012).}
  \item{log.p}{\code{\link{logical}}; if \code{TRUE}, probabilities
    \eqn{p} are given as \eqn{\log p}{log(p)}.}
  \item{p}{probabilities or log-probabilities if \code{log.p} is true.}
  \item{method}{string for the method to compute the quantile function
    of \eqn{K}. Currently, one of
    \describe{
      \item{"default"}{default method. Currently chooses
	\code{method="monoH.FC"} with \code{u.grid = 0:128/128}. This
	is fast but not too accurate (see example).}
      \item{"simple"}{straightforward root finding based on
	\code{\link{uniroot}}.}
      \item{"sort"}{root finding based on \code{\link{uniroot}} but
	first sorting \code{u}.}
      \item{"discrete"}{first, \eqn{K} is evaluated at the given grid
	points \code{u.grid} (which should contain 0 and 1). Based on these
	probabilities, quantiles are computed with \code{\link{findInterval}}.}
      \item{"monoH.FC"}{first, \eqn{K} is evaluated at the given grid points
	\code{u.grid}.  A monotone spline is then used to approximate \eqn{K}.
	Based on this approximation, quantiles are computed with
	\code{\link{uniroot}}.}
    }
  }
  \item{u.grid}{(for \code{method="discrete"}:) The grid on which
    \eqn{K} is evaluated, a \code{\link{numeric}} vector.}
  \item{\dots}{additional arguments passed to \code{\link{uniroot}} (for
    \code{method="default"}, \code{method="simple"}, \code{method="sort"}, and
    \code{method="monoH.FC"}) or \code{\link{findInterval}} (for
    \code{method="discrete"}), notably \code{tol}
    (\code{\link{uniroot}}) for increased accuracy.}
  \item{n}{sample size for \code{rK}.}
}
\details{
  For a completely monotone Archimedean generator \eqn{\psi}{psi},
  \deqn{K(u)=\sum_{k=0}^{d-1}
    \frac{\psi^{(k)}(\psi^{-1}(u))}{k!} (-\psi^{-1}(u))^k,\ u\in[0,1];}{%
    K(u)=sum(k=0,...,d-1) psi^{(k)}(psi^{-1}(u))/k! (-psi^{-1}(u))^k, u
    in [0,1];}
  see Barbe et al. (1996).  The corresponding density is
  \deqn{\frac{(-1)^d\psi^{(d)}(\psi^{-1}(u))}{(d-1)!}
    (-(\psi^{-1})'(u))(\psi^{-1}(u))^{d-1}}{%
    (-1)^d
    psi^{(d)}(psi^{-1}(u))(-(psi^{-1})'(u))(psi^{-1}(u))^{d-1}/(d-1)!}

  The empirical Kendall distribution function is computed as in Genest,
  G. \enc{Nešlehová}{Neslehova}, Ziegel (2011).
}
\value{
  The empirical Kendall distribution function, density, distribution
  function, quantile function and random number generator.
}
%\author{Marius Hofert (and Martin Maechler)}
\note{
  Currently, the \code{"default"} method of \code{qK()} is fast but
  not very accurate, see the \sQuote{Examples} for more accuracy (with
  more CPU effort).
}
\references{
  Barbe, P., Genest, C., Ghoudi, K., and R\enc{é}{e}millard, B. (1996),
  On Kendall's Process,
  \emph{Journal of Multivariate Analysis} \bold{58}, 197--229.

  Hofert, M., \enc{Mächler}{Maechler}, M., and McNeil, A. J. (2012).
  Likelihood inference for Archimedean copulas in high dimensions
  under known margins. \emph{Journal of Multivariate Analysis}
  \bold{110}, 133--150.

  Genest, C., G. \enc{Nešlehová}{Neslehova}, J., and Ziegel, J. (2011).
  Inference in multivariate Archimedean copula models. \emph{TEST}
  \bold{20}, 223--256.
}
\seealso{
  \code{\link{htrafo}} or \code{\link{emde}} (where \code{\link{K}} is used);
  \code{\link{splinefun}(*, "monoHC")} for that method.
}
\examples{
tau <- 0.5
(theta <- copGumbel@iTau(tau)) # 2
d <- 20
(cop <- onacopulaL("Gumbel", list(theta,1:d)))

## Basic check of the empirical Kendall distribution function
set.seed(271)
n <- 1000
U <- rCopula(n, copula = cop)
X <- qnorm(U)
K.sample <- pCopula(U, copula = cop)
u <- seq(0, 1, length.out = 256)
edfK <- ecdf(K.sample)
plot(u, edfK(u), type = "l", ylim = 0:1,
     xlab = quote(italic(u)), ylab = quote(K[n](italic(u)))) # simulated
K.n <- Kn(u, x = X)
lines(u, K.n, col = "royalblue3") # Kn
## Difference at 0
edfK(0) # edf of K at 0
K.n[1] # K_n(0); this is > 0 since K.n is the edf of a discrete distribution
## => therefore, Kn(K.sample, x=X) is not uniform
plot(Kn(K.sample, x = X), ylim = 0:1)
## Note: Kn(0) -> 0 for n -> Inf

## Compute Kendall distribution function
u <- seq(0,1, length.out = 255)
Ku    <- pK(u, copula = cop@copula, d = d) # exact
Ku.MC <- pK(u, copula = cop@copula, d = d, n.MC = 1000) # via Monte Carlo
stopifnot(all.equal(log(Ku),
		    pK(u, copula = cop@copula, d = d, log.p=TRUE)))# rel.err 3.2e-16

## Build sample from K
set.seed(1)
n <- 200
W <- rK(n, copula = cop)

## Plot empirical distribution function based on W
## and the corresponding theoretical Kendall distribution function
## (exact and via Monte Carlo)
plot(ecdf(W), col = "blue", xlim = 0:1, verticals=TRUE,
     main = quote("Empirical"~ F[n](C(U)) ~
                     "and its Kendall distribution" ~ K(u)),
     do.points = FALSE, asp = 1)
abline(0,1, lty = 2); abline(h = 0:1, v = 0:1, lty = 3, col = "gray")
lines(u, Ku.MC, col = "red") # not quite monotone
lines(u, Ku, col = "black")  # strictly  monotone:
stopifnot(diff(Ku) >= 0)
legend(.25, .75, expression(F[n], K[MC](u), K(u)),
       col=c("blue" , "red", "black"), lty = 1, lwd = 1.5, bty = "n")

if(require("Rmpfr")) { # pK() now also works with high precision numbers:
 uM <- mpfr(0:255, 99)/256
 if(FALSE) {
   # not yet, now fails in  polyG() :
   KuM <- pK(uM, copula = cop@copula, d = d)
  ##  debug(copula:::.pK)
  debug(copula:::polyG)
 }
}# if( Rmpfr )


## Testing qK
pexpr <- quote( 0:63/63 );  p <- eval(pexpr)
d <- 10
cop <- onacopulaL("Gumbel", list(theta = 2, 1:d))
system.time(qK0 <- qK(p, copula = cop@copula, d = d)) # "default" - fast
\donttest{% needs a couple of valuable CRAN seconds :-)

system.time(qK1  <- qK(p, copula= cop@copula, d=d, method = "simple"))
system.time(qK1. <- qK(p, copula= cop@copula, d=d, method = "simple", tol = 1e-12))
system.time(qK2  <- qK(p, copula= cop@copula, d=d, method = "sort"))
system.time(qK2. <- qK(p, copula= cop@copula, d=d, method = "sort",   tol = 1e-12))
system.time(qK3  <- qK(p, copula= cop@copula, d=d, method = "discrete", u.grid = 0:1e4/1e4))
system.time(qK4  <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = 0:5e2/5e2))
system.time(qK4. <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = 0:5e2/5e2, tol = 1e-12))
system.time(qK5  <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = 0:5e3/5e3))
system.time(qK5. <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = 0:5e3/5e3, tol = 1e-12))
system.time(qK6  <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = (0:5e3/5e3)^2))
system.time(qK6. <- qK(p, copula= cop@copula, d=d, method = "monoH.FC",
                       u.grid = (0:5e3/5e3)^2, tol = 1e-12))

## Visually they all coincide :
cols <- adjustcolor(c("gray50", "gray80", "light blue",
                      "royal blue", "purple3", "purple4", "purple"), 0.6)
matplot(p, cbind(qK0, qK1, qK2, qK3, qK4, qK5, qK6), type = "l", lwd = 2*7:1, lty = 1:7, col = cols,
        xlab = bquote(p == .(pexpr)), ylab = quote({K^{-1}}(u)),
        main = "qK(p, method = *)")
legend("topleft", col = cols, lwd = 2*7:1, lty = 1:7, bty = "n", inset = .03,
       legend= paste0("method= ",
             sQuote(c("default", "simple", "sort",
                      "discrete(1e4)", "monoH.FC(500)", "monoH.FC(5e3)", "monoH.FC(*^2)"))))

## See they *are* inverses  (but only approximately!):
eqInv <- function(qK) all.equal(p, pK(qK, cop@copula, d=d), tol=0)

eqInv(qK0 ) # "default"	       0.03  worst
eqInv(qK1 ) # "simple"	       0.0011 - best
eqInv(qK1.) # "simple", e-12   0.00000 (8.73 e-13) !
eqInv(qK2 ) # "sort"	       0.0013 (close)
eqInv(qK2.) # "sort", e-12     0.00000 (7.32 e-12)
eqInv(qK3 ) # "discrete"       0.0026
eqInv(qK4 ) # "monoH.FC(500)"  0.0095
eqInv(qK4.) # "m.H.FC(5c)e-12" 0.00963
eqInv(qK5 ) # "monoH.FC(5e3)"  0.001148
eqInv(qK5.) # "m.H.FC(5k)e-12" 0.000989
eqInv(qK6 ) # "monoH.FC(*^2)"  0.001111
eqInv(qK6.) # "m.H.FC(*^2)e-12"0.00000 (1.190 e-09)

## and ensure the differences are not too large
stopifnot(
 all.equal(qK0, qK1, tol = 1e-2) # !
 ,
 all.equal(qK1, qK2, tol = 1e-4)
 ,
 all.equal(qK2, qK3, tol = 1e-3)
 ,
 all.equal(qK3, qK4, tol = 1e-3)
 ,
 all.equal(qK4, qK0, tol = 1e-2) # !
)

}%dont
stopifnot(all.equal(p, pK(qK0, cop@copula, d=d), tol = 0.04))
}
\keyword{distribution}
