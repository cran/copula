\name{Copula}
\alias{Copula}
\alias{dcopula}
\alias{pcopula}
\alias{rcopula}

\alias{dcopula,claytonCopula-method}
\alias{dcopula,frankCopula-method}
\alias{dcopula,gumbelCopula-method}
\alias{dcopula,normalCopula-method}
\alias{dcopula,tCopula-method}
\alias{dcopula,galambosCopula-method}
\alias{dcopula,amhCopula-method}
\alias{dcopula,huslerReissCopula-method}
\alias{dcopula,tawnCopula-method}
\alias{dcopula,tevCopula-method}
\alias{dcopula,plackettCopula-method}
\alias{dcopula,fgmCopula-method}
\alias{dcopula,indepCopula-method}
\alias{dcopula,asymCopula-method}
\alias{dcopula,asymExplicitCopula-method}

\alias{pcopula,claytonCopula-method}
\alias{pcopula,frankCopula-method}
\alias{pcopula,gumbelCopula-method}
\alias{pcopula,normalCopula-method}
\alias{pcopula,tCopula-method}
\alias{pcopula,galambosCopula-method}
\alias{pcopula,amhCopula-method}
\alias{pcopula,huslerReissCopula-method}
\alias{pcopula,tawnCopula-method}
\alias{pcopula,tevCopula-method}
\alias{pcopula,plackettCopula-method}
\alias{pcopula,fgmCopula-method}
\alias{pcopula,indepCopula-method}
\alias{pcopula,asymCopula-method}
\alias{pcopula,asymExplicitCopula-method}


\alias{rcopula,claytonCopula-method}
\alias{rcopula,frankCopula-method}
\alias{rcopula,gumbelCopula-method}
\alias{rcopula,normalCopula-method}
\alias{rcopula,tCopula-method}
\alias{rcopula,galambosCopula-method}
\alias{rcopula,amhCopula-method}
\alias{rcopula,huslerReissCopula-method}
\alias{rcopula,plackettCopula-method}
\alias{rcopula,fgmCopula-method}
\alias{rcopula,evCopula-method}
\alias{rcopula,indepCopula-method}
\alias{rcopula,asymCopula-method}
\alias{rcopula,asymExplicitCopula-method}

\title{Copula distribution functions}
\description{
  Density, distribution function, and random generation for a
  \code{copula} object.
}
\usage{
dcopula(copula, u)
pcopula(copula, u)
rcopula(copula, n)
}
\arguments{
  \item{copula}{ a \code{"copula"} object. }
  \item{u}{ a vector of the copula dimension or a matrix with number of
    rows being the copula dimension, giving the coordinates of the
    points where the density or distribution function needs to be evaluated.
  }
  \item{n}{ number of observations to be generated. }
}
\details{
  The density function of an Archimedean copula was obtained by
  differentiating the distribution function symbolically using
  Mathematica and then processed by \code{deriv} to give
  algorithmic expressions. The maximum dimension implemented is 10 for
  Clayton and Gumbel, and 6 for Frank.
  
  The distribution function of a \code{t} copula uses \code{pmvt} from
  package \code{mvtnorm}. The density function of a \code{t} copula uses the
  \code{dmst} from package \code{sn}.
  
  The random number generator for an Archimedean copula uses the
  conditional approach for the bivariate case and the Marshal-Olkin (1988)
  approach for dimension greater than 2.
}
\value{
  \code{"dcopula"} gives the density, \code{"pcopula"} gives the
  distribution function, and \code{"rcopula"} generates random variates.
}
\references{
  E.W. Frees and E.A. Valdez (1998). Understanding relationships using
  copulas.  \emph{North American Actuarial Journal}, 2:1--25.
  
  C. Genest and A.-C. Favre (2007). Everything you always wanted to know
  about copula modeling but were afraid to ask. \emph{Journal of
    Hydrologic Engineering}, 12:347--368.
  
  H. Joe (1997). \emph{Multivariate Models and Dependence Concepts}.
  Chapman and Hall, London.

  A.W. Marshal and I. Olkin (1988). Families of multivariate
  distributions. \emph{Journal of the American Statistical Association},
  83:834--841.
      
  R.B. Nelsen (2006). \emph{An introduction to Copulas}. Springer, New York.
}
\seealso{
  \code{\link{copula-class}},
  \code{\link{ellipCopula}},
  \code{\link{archmCopula}},
  \code{\link{fgmCopula}}. }
\examples{
norm.cop <- normalCopula(0.5)
norm.cop
x <- rcopula(norm.cop, 100)
plot(x)
dcopula(norm.cop, x)
pcopula(norm.cop, x)
persp(norm.cop, dcopula)
contour(norm.cop, pcopula)
## a 3-dimensional normal copula
u <- rcopula(normalCopula(0.5, dim = 3), 1000)
## scatterplot3d(u)
## a 3-dimensional clayton copula
v <- rcopula(claytonCopula(2, dim = 3), 1000)
## scatterplot3d(v)
}
\keyword{distribution}% at least one, from doc/KEYWORDS
\keyword{multivariate}% __ONLY ONE__ keyword per line
