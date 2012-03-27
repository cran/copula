\name{nacPairthetas}
\alias{nacPairthetas}
\title{Pairwise thetas of Nested Archimedean Copula}
\description{
  Return a \eqn{d * d} matrix of pairwise thetas for a nested Archimedean
  copula (\code{\linkS4class{nacopula}}) of dimension \eqn{d}.
}
\usage{
nacPairthetas(x)
}
\arguments{
  \item{x}{an (outer) nacopula (with thetas sets).}
}
% \details{
% %%  ~~ If necessary, more details than the description above ~~
% }
\value{
  a (\eqn{d \times d}{d x d}) matrix of thetas, say \code{T}, where
  \code{T[j,k]} = theta of the bivariate Archimedean copula
  \eqn{C(U_j,U_k)}.
}
\author{Martin Maechler}
% \references{
% }
\seealso{
  the class \code{\linkS4class{nacopula}} (with its \code{\link{dim}}
  method).
}
\examples{
## test with
options(width=97)

(mm <- rnacModel("Gumbel", d=15, pr.comp = 0.25, order="random"))
stopifnot(isSymmetric(PT <- nacPairthetas(mm)))
round(PT, 2)

## do this several times:
m1 <- rnacModel("Gumbel", d=15, pr.comp = 1/8, order="seq")
stopifnot(isSymmetric(PT <- nacPairthetas(m1)))
m1; PT

%%___ FIXME __: This shows that rnacModel() should  choose  'nkids' less aggressively..
m100 <- rnacModel("Gumbel", d= 100, pr.comp = 1/16, order="seq")
system.time(PT <- nacPairthetas(m100))# how slow {non-optimal algorithm} ??
##-- very fast, still!
stopifnot(isSymmetric(PT))
m100

## image(PT)# not ok -- want one color per theta
nt <- length(th0 <- unique(PT[!is.na(PT)]))
th1 <- c(th0[1]/2, th0, 1.25*th0[nt])
ths <- (th1[-1]+th1[-(nt+2)])/2
image(log(PT), breaks = ths, col = heat.colors(nt))

## Nicer and easier:
require("Matrix")
image(as(log(PT),"Matrix"), main = "log( nacPairthetas( m100 ))",
      useAbs=FALSE, useRaster=TRUE, border=NA)
}
\keyword{utilities}
