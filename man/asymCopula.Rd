\name{khoudrajiCopula}
\title{Construction of copulas using Khoudraji's device}
\alias{khoudrajiCopula}
% Deprecated:
\alias{asymCopula}
\alias{asymExplicitCopula}
\description{
  Creates an object representing a copula constructed using
  \emph{Khoudraji's device} (Khoudraji, 1995). The resulting \R object
  is either of class \code{"\linkS4class{khoudrajiBivCopula}"},
  \code{"\linkS4class{khoudrajiExplicitCopula}"}
  or \code{"\linkS4class{khoudrajiCopula}"}.

  In the bivariate case, given two copulas \eqn{C_1} and \eqn{C_2},
  Khoudraji's device consists of defining a copula whose c.d.f. is given
  by:

  \deqn{C_1(u_1^{1 - a_1}, u_2^{1 - a_2}) C_2(u_1^{a_1}, u_2^{a_2})}{
    C_1(u_1^{1 - a_1}, u_2^{1 - a_2}) C_2(u_1^{a_1}, u_2^{a_2})}

  where \eqn{a_1} and \eqn{a_2} are \emph{shape parameters} in [0,1].

  The construction principle (see also Genest et al. 1998) is a special
  case of that considered in Liebscher (2008).
}
\usage{
khoudrajiCopula(copula1 = indepCopula(), copula2 = indepCopula(),
                shapes = rep(NA_real_, dim(copula1)))
}
\arguments{
  \item{copula1, copula2}{each a \code{\linkS4class{copula}}, of the
    same dimension \eqn{d}.}
  \item{shapes}{\code{\link{numeric}} vector of length \eqn{d},
    with values in \eqn{[0,1]}.}
}
\details{

  If the argument copulas are bivariate, an object of class
  \code{"\linkS4class{khoudrajiBivCopula}"} will be constructed.  If
  they are exchangeable and \eqn{d}-dimensional with \eqn{d > 2}, and if
  they have explicit p.d.f. and c.d.f. expressions, an object of class
  \code{"\linkS4class{khoudrajiExplicitCopula}"} will be
  constructed. For the latter two classes, density evaluation is implemented,
  and fitting and goodness-of-fit testing can be attempted. If \eqn{d >
  2} but one of the argument copulas does not have explicit p.d.f. and
  c.d.f. expressions, or is not exchangeable, an object of class
  \code{"\linkS4class{khoudrajiCopula}"} will be constructed, for which
  density evaluation is not possible.
}
%\author{Jun Yan and Ivan Kojadinovic}
\value{
  A new object of class \code{"\linkS4class{khoudrajiBivCopula}"} in
  dimension two or of class
  \code{"\linkS4class{khoudrajiExplicitCopula}"} or
  \code{"\linkS4class{khoudrajiCopula}"} when \eqn{d > 2}.
}
\references{
  Genest, C., Ghoudi, K., and Rivest, L.-P. (1998), Discussion of
  "Understanding relationships using copulas", by Frees, E., and
  Valdez, E., \emph{North American Actuarial Journal} \bold{3},
  143--149.

  Khoudraji, A. (1995), Contributions à l'\enc{é}{e}tude des copules et
  \enc{à}{a}la mod\enc{é}{e}lisation des valeurs extr\enc{ê}{e}mes
  bivari\enc{é}{e}es, \emph{PhD thesis, Université Laval}, Québec,
  Canada.

  Liebscher, E. (2008), Construction of asymmetric multivariate copulas,
  \emph{Journal of Multivariate Analysis} \bold{99}, 2234--2250.
}
%% \seealso{
%% }
\examples{
## A bivariate Khoudraji-Clayton copula
kc <- khoudrajiCopula(copula2 = claytonCopula(6),
                      shapes = c(0.4, 0.95))
class(kc)

kc ## FIXME: too long 'fullname'
contour(kc, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

## A Khoudraji-Clayton copula with second shape parameter fixed
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))

## A "nested" Khoudraji bivariate copula
kgkcf <- khoudrajiCopula(copula1 = gumbelCopula(3),
                         copula2 = kcf,
                         shapes = c(0.7, 0.25))
contour(kgkcf, dCopula, nlevels = 20, main = "dCopula(<khoudrajiBivCopula>)")

n <- 300
u <- rCopula(n, kc)
plot(u)

\dontrun{
## Fitting can be difficult
## Starting values are required for all parameters
fitCopula(khoudrajiCopula(copula2 = claytonCopula()),
          start = c(1.1, 0.5, 0.5), data = pobs(u),
          optim.method = "Nelder-Mead")

## Second shape parameter fixed to 0.95
kcf2 <- khoudrajiCopula(copula2 = claytonCopula(),
                        shapes = fixParam(c(NA_real_, 0.95), c(FALSE, TRUE)))
fitCopula(kcf2,
          start = c(1.1, 0.5), data = pobs(u),
          optim.method = "Nelder-Mead")

## With a different optimization method
fitCopula(kcf2,
          start = c(1.1, 0.5), data = pobs(u),
          optim.method = "BFGS")

## GOF example
optim.method <- "Nelder-Mead" #try "BFGS" as well
gofCopula(kcf2, x = u, start = c(1.1, 0.5), optim.method = optim.method)
gofCopula(kcf2, x = u, start = c(1.1, 0.5), optim.method = optim.method,
          sim = "mult")
## The goodness-of-fit tests should hold their level
## but this would need to be tested

## Another example under the alternative
u <- rCopula(n, gumbelCopula(4))
gofCopula(kcf2, x = u, start = c(1.1, 0.5), optim.method = optim.method)
gofCopula(kcf2, x = u, start = c(1.1, 0.5), optim.method = optim.method,
          sim = "mult")

## Higher-dimensional constructions

## A three dimensional Khoudraji-Clayton copula
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3
class(kcd3)

n <- 1000
u <- rCopula(n, kcd3)
v <- matrix(runif(15), 5, 3)
splom2(u)
dCopula(v, kcd3)

## A four dimensional Khoudraji-Normal copula
knd4 <- khoudrajiCopula(copula1 = indepCopula(dim=4),
                        copula2 = normalCopula(.9, dim=4),
                        shapes = c(0.4, 0.95, 0.95, 0.95))
knd4
class(knd4)
u <- rCopula(n, knd4)
splom2(u)
## dCopula(v, knd4) ## not implemented
}
}

\keyword{distribution}
\keyword{multivariate}
