\name{Copula}
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
\alias{dCopula,matrix,indepCopula-method}
\alias{dCopula,matrix,asymCopula-method}
\alias{dCopula,matrix,asymExplicitCopula-method}
\alias{dCopula,matrix,nacopula-method}
%
\alias{dCopula,numeric,claytonCopula-method}
\alias{dCopula,numeric,frankCopula-method}
\alias{dCopula,numeric,gumbelCopula-method}
\alias{dCopula,numeric,normalCopula-method}
\alias{dCopula,numeric,tCopula-method}
\alias{dCopula,numeric,galambosCopula-method}
\alias{dCopula,numeric,amhCopula-method}
\alias{dCopula,numeric,joeCopula-method}
\alias{dCopula,numeric,huslerReissCopula-method}
\alias{dCopula,numeric,tawnCopula-method}
\alias{dCopula,numeric,tevCopula-method}
\alias{dCopula,numeric,plackettCopula-method}
\alias{dCopula,numeric,fgmCopula-method}
\alias{dCopula,numeric,indepCopula-method}
\alias{dCopula,numeric,asymCopula-method}
\alias{dCopula,numeric,asymExplicitCopula-method}
\alias{dCopula,numeric,nacopula-method}
% pCopula
\alias{pCopula,matrix,claytonCopula-method}
\alias{pCopula,matrix,frankCopula-method}
\alias{pCopula,matrix,gumbelCopula-method}
\alias{pCopula,matrix,normalCopula-method}
\alias{pCopula,matrix,tCopula-method}
\alias{pCopula,matrix,galambosCopula-method}
\alias{pCopula,matrix,amhCopula-method}
\alias{pCopula,matrix,joeCopula-method}
\alias{pCopula,matrix,huslerReissCopula-method}
\alias{pCopula,matrix,tawnCopula-method}
\alias{pCopula,matrix,tevCopula-method}
\alias{pCopula,matrix,plackettCopula-method}
\alias{pCopula,matrix,fgmCopula-method}
\alias{pCopula,matrix,indepCopula-method}
\alias{pCopula,matrix,asymCopula-method}
\alias{pCopula,matrix,asymExplicitCopula-method}
\alias{pCopula,matrix,nacopula-method}
%
\alias{pCopula,numeric,claytonCopula-method}
\alias{pCopula,numeric,frankCopula-method}
\alias{pCopula,numeric,gumbelCopula-method}
\alias{pCopula,numeric,normalCopula-method}
\alias{pCopula,numeric,tCopula-method}
\alias{pCopula,numeric,galambosCopula-method}
\alias{pCopula,numeric,amhCopula-method}
\alias{pCopula,numeric,joeCopula-method}
\alias{pCopula,numeric,huslerReissCopula-method}
\alias{pCopula,numeric,tawnCopula-method}
\alias{pCopula,numeric,tevCopula-method}
\alias{pCopula,numeric,plackettCopula-method}
\alias{pCopula,numeric,fgmCopula-method}
\alias{pCopula,numeric,indepCopula-method}
\alias{pCopula,numeric,asymCopula-method}
\alias{pCopula,numeric,asymExplicitCopula-method}
\alias{pCopula,numeric,nacopula-method}
% rCopula
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
\alias{rCopula,numeric,evCopula-method}
\alias{rCopula,numeric,indepCopula-method}
\alias{rCopula,numeric,asymCopula-method}
\alias{rCopula,numeric,asymExplicitCopula-method}
\alias{rCopula,numeric,nacopula-method}
% main
\title{Density, Evaluation, and Random Number Generation for Copula Functions}
\description{
  Density, distribution function, and random generation for a
  \code{copula} object.
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
    function needs to be evaluated.
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
u <- rCopula(100, norm.cop)
plot(u)
dCopula(u, norm.cop)
pCopula(u, norm.cop)
persp(norm.cop, dCopula)
contour(norm.cop, pCopula)

## a 3-dimensional normal copula
u <- rCopula(1000, normalCopula(0.5, dim = 3))
if(require("scatterplot3d"))
  scatterplot3d(u)

## a 3-dimensional clayton copula
cl3 <- claytonCopula(2, dim = 3)
v <- rCopula(1000, cl3)
pairs(v)
if(require("scatterplot3d"))
  scatterplot3d(v)

## Compare with the "nacopula" version :
fu1 <- dCopula(v, cl3)
fu2 <- copClayton@dacopula(v, theta = 2)
Fu1 <- pCopula(v, cl3)
Fu2 <- pnacopula(onacopula("Clayton", C(2.0, 1:3)), v)
## The density and cumulative values are the same:
stopifnot(all.equal(fu1, fu2, tol= 1e-14),
          all.equal(Fu1, Fu2, tol= 1e-15))
}
\keyword{distribution}
\keyword{multivariate}
