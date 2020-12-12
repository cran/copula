\name{Copula}
\title{Density, Evaluation, and Random Number Generation for Copula Functions}
\alias{Copula}
\alias{dCopula}
\alias{pCopula}
\alias{rCopula}
\alias{dcopula}%deprecated
\alias{pcopula}%deprecated
\alias{rcopula}%deprecated
% dCopula
\alias{dCopula,matrix,claytonCopula-method}
\alias{dCopula,matrix,frankCopula-method}
\alias{dCopula,matrix,gumbelCopula-method}
\alias{dCopula,matrix,normalCopula-method}
\alias{dCopula,matrix,tCopula-method}
\alias{dCopula,matrix,galambosCopula-method}
\alias{dCopula,matrix,amhCopula-method}
\alias{dCopula,matrix,joeCopula-method}
\alias{dCopula,matrix,huslerReissCopula-method}
\alias{dCopula,matrix,tawnCopula-method}
\alias{dCopula,matrix,tevCopula-method}
\alias{dCopula,matrix,plackettCopula-method}
\alias{dCopula,matrix,fgmCopula-method}
\alias{dCopula,matrix,moCopula-method}
\alias{dCopula,matrix,indepCopula-method}
\alias{dCopula,matrix,fhCopula-method}
\alias{dCopula,matrix,empCopula-method}
\alias{dCopula,matrix,khoudrajiBivCopula-method}
\alias{dCopula,matrix,khoudrajiExplicitCopula-method}
\alias{dCopula,matrix,rotCopula-method}
\alias{dCopula,matrix,rotExplicitCopula-method}
\alias{dCopula,matrix,mixCopula-method}
%\alias{dCopula,matrix,nacopula-method}-> ./dnacopula.Rd
%
% pCopula ====
\alias{pCopula,matrix,claytonCopula-method}
\alias{pCopula,matrix,frankCopula-method}
\alias{pCopula,matrix,gumbelCopula-method}
\alias{pCopula,matrix,galambosCopula-method}
\alias{pCopula,matrix,amhCopula-method}
\alias{pCopula,matrix,joeCopula-method}
\alias{pCopula,matrix,huslerReissCopula-method}
%\alias{pCopula,matrix,normalCopula-method}% -> ./ellipCopula.Rd
%\alias{pCopula,matrix,tCopula-method}     % ->     (ditto)
\alias{pCopula,matrix,tawnCopula-method}
\alias{pCopula,matrix,tevCopula-method}
\alias{pCopula,matrix,plackettCopula-method}
\alias{pCopula,matrix,fgmCopula-method}
\alias{pCopula,matrix,moCopula-method}
\alias{pCopula,matrix,indepCopula-method}
\alias{pCopula,matrix,lowfhCopula-method}
\alias{pCopula,matrix,upfhCopula-method}
\alias{pCopula,matrix,empCopula-method}
\alias{pCopula,matrix,khoudrajiCopula-method}
\alias{pCopula,matrix,rotCopula-method}
\alias{pCopula,matrix,rotExplicitCopula-method}
%\alias{pCopula,matrix,nacopula-method}-> ./pnacopula.Rd
\alias{pCopula,matrix,mixCopula-method}
%
% rCopula ====
\alias{rCopula,numeric,claytonCopula-method}
\alias{rCopula,numeric,frankCopula-method}
\alias{rCopula,numeric,gumbelCopula-method}
\alias{rCopula,numeric,normalCopula-method}
\alias{rCopula,numeric,tCopula-method}
\alias{rCopula,numeric,galambosCopula-method}
\alias{rCopula,numeric,amhCopula-method}
\alias{rCopula,numeric,joeCopula-method}
\alias{rCopula,numeric,huslerReissCopula-method}
\alias{rCopula,numeric,plackettCopula-method}
\alias{rCopula,numeric,fgmCopula-method}
\alias{rCopula,numeric,moCopula-method}
\alias{rCopula,numeric,evCopula-method}
\alias{rCopula,numeric,indepCopula-method}
\alias{rCopula,numeric,lowfhCopula-method}
\alias{rCopula,numeric,upfhCopula-method}
\alias{rCopula,numeric,empCopula-method}
\alias{rCopula,numeric,khoudrajiCopula-method}
\alias{rCopula,numeric,rotCopula-method}
\alias{rCopula,numeric,mixCopula-method}
\alias{rCopula,numeric,nacopula-method}
% main
\description{
  Density (\code{dCopula}), distribution function (\code{pCopula}), and
  random generation (\code{rCopula}) for a \code{copula} object.
}
\usage{
dCopula(u, copula, log=FALSE, \dots)
pCopula(u, copula, \dots)
rCopula(n, copula, \dots)
}
\arguments{
  \item{copula}{an \R object of class \code{"\linkS4class{Copula}"},
    (i.e., \code{"\linkS4class{copula}"} or  \code{"\linkS4class{nacopula}"}).}
  \item{u}{a vector of the copula dimension \eqn{d} or a matrix with \eqn{d}
    columns, giving the points where the density or distribution
    function needs to be evaluated.  Note that in all cases, values
    outside of the cube \eqn{[0,1]^d} are treated equivalently
    to those on the cube boundary.  So, e.g., the density is zero.
  }
  \item{log}{logical indicating if the \eqn{\log(f(\cdot))}{log(f(.))}
    should be returned instead of \eqn{f(\cdot)}{f(.)}.}
  \item{n}{(for \code{rCopula()}:) number of observations to be generated.}
  \item{\dots}{further optional arguments for some methods, e.g.,
    \code{method}.}
}
\details{
  The density (\code{dCopula}) and distribution function
  (\code{pCopula}) methods for Archimedean copulas now use the
  corresponding function slots of the Archimedean copula objects, such
  as \code{\link{copClayton}}, \code{\link{copGumbel}}, etc.

  If an \eqn{u_j}{u[j]} is outside \eqn{(0,1)} we declare the density to be
  zero, and this is true even when another \eqn{u_k, k \ne j}{u[k], k != j}
  is \code{\link{NA}} or \code{\link{NaN}}; see also the
  \dQuote{outside} example.

  The distribution function of a \code{t} copula uses \code{pmvt} from
  package \pkg{mvtnorm}; similarly, the density (\code{dCopula}) calls
  \code{\link[mvtnorm]{dmvt}} from \pkg{mvtnorm}.
  The \code{normalCopula} methods use \code{dmvnorm} and \code{pmvnorm}
  from the same package.

  The random number generator for an Archimedean copula uses the
  conditional approach for the bivariate case and the Marshall-Olkin (1988)
  approach for dimension greater than 2.
}
\value{
  \code{dCopula()} gives the density, \code{pCopula()} gives the
  distribution function, and \code{rCopula()} generates random variates.
}
%\author{Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan}
\references{
  Frees, E. W.  and Valdez, E. A. (1998). Understanding relationships using
  copulas.
  \emph{North American Actuarial Journal} \bold{2}, 1--25.

  Genest, C. and Favre, A.-C. (2007).
  Everything you always wanted to know about copula modeling but were
  afraid to ask.
  \emph{Journal of Hydrologic Engineering} \bold{12}, 347--368.

  Joe, H. (1997). \emph{Multivariate Models and Dependence Concepts}.
  Chapman and Hall, London.

  Marshall, A. W. and Olkin, I. (1988)
  Families of multivariate distributions.
  \emph{Journal of the American Statistical Association} \bold{83}, 834--841.

  Nelsen, R. B. (2006)
  \emph{An introduction to Copulas}. Springer, New York.
}
\seealso{
  the \code{\linkS4class{copula}} and \code{\linkS4class{acopula}}
  classes, the acopula families, \code{\link{acopula-families}}.
  Constructor functions such as
  \code{\link{ellipCopula}}, \code{\link{archmCopula}},
  \code{\link{fgmCopula}}.
}
\examples{
norm.cop <- normalCopula(0.5)
norm.cop
## one d-vector =^= 1-row matrix, works too :
dCopula(c(0.5, 0.5), norm.cop)
pCopula(c(0.5, 0.5), norm.cop)

u <- rCopula(100, norm.cop)
plot(u)
dCopula(u, norm.cop)
pCopula(u, norm.cop)
persp  (norm.cop, dCopula)
contour(norm.cop, pCopula)

## a 3-dimensional normal copula
u <- rCopula(1000, normalCopula(0.5, dim = 3))
if(require(scatterplot3d))
  scatterplot3d(u)

## a 3-dimensional clayton copula
cl3 <- claytonCopula(2, dim = 3)
v <- rCopula(1000, cl3)
pairs(v)
if(require(scatterplot3d))
  scatterplot3d(v)

## Compare with the "nacopula" version :
fu1 <- dCopula(v, cl3)
fu2 <- copClayton@dacopula(v, theta = 2)
Fu1 <- pCopula(v, cl3)
Fu2 <- pCopula(v, onacopula("Clayton", C(2.0, 1:3)))
## The density and cumulative values are the same:
stopifnot(all.equal(fu1, fu2, tolerance= 1e-14),
          all.equal(Fu1, Fu2, tolerance= 1e-15))

## NA and "outside" u[]
u <- v[1:12,]
## replace some by values outside (0,1) and some by NA/NaN
u[1, 2:3] <- c(1.5, NaN); u[2, 1] <- 2; u[3, 1:2] <- c(NA, -1)
u[cbind(4:9, 1:3)] <- c(NA, NaN)
f <- dCopula(u, cl3)
cbind(u, f) # note: f(.) == 0 at [1] and [3] inspite of NaN/NA
stopifnot(f[1:3] == 0, is.na(f[4:9]), 0 < f[10:12])
}
\keyword{distribution}
\keyword{multivariate}
