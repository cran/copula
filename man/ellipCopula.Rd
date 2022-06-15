\name{ellipCopula}
\title{Construction of Elliptical Copula Class Objects}
\alias{ellipCopula}
\alias{normalCopula}
\alias{tCopula}
\alias{dispstrToep}
\alias{pCopula,matrix,normalCopula-method}
\alias{pCopula,matrix,tCopula-method}
\description{
  Creating elliptical copula objects with corresponding dimension and
  parameters, including the dispersion structure \eqn{P} (pronounced \dQuote{Rho}).
}
\usage{
ellipCopula (family, param, dim = 2, dispstr = "ex", df = 4, ...)

normalCopula(param, dim = 2, dispstr = "ex")

     tCopula(param, dim = 2, dispstr = "ex", df = 4,
             df.fixed = FALSE, df.min = 0.01)

dispstrToep(perm = NULL, check = TRUE)

\S4method{pCopula}{matrix,normalCopula}(u, copula, algorithm=NULL, keepAttr=FALSE, ...)
% --> pnormalCopula() in ../R/normalCopula.R
\S4method{pCopula}{matrix,tCopula}(u, copula, algorithm=NULL, keepAttr=FALSE, ...)
% --> ptCopula() in ../R/tCopula.R
}
\arguments{
  \item{family}{a \code{\link{character}} string specifying the family of an
    elliptical copula.  Must be \code{"normal"} (the default) or \code{"t"}.}
  \item{param}{a \code{\link{numeric}} vector specifying the parameter values;
    \code{P2p()} accesses this vector, whereas
    \code{\link{p2P}()} and \code{\link{getSigma}()} provide
    the corresponding \dQuote{P} matrix, see below.}
  \item{dim}{ the dimension of the copula. }
  \item{dispstr}{a string specifying the \dQuote{dispersion structure},
    i.e., type of the symmetric positive definite matrix characterizing the
    elliptical copula.  Currently available structures are \code{"ex"} for
    \bold{ex}changeable, \code{"ar1"} for \eqn{AR(1)}, \code{"toep"} for
    Toeplitz (\code{\link{toeplitz}}), and \code{"un"} for
    \bold{un}structured.

    The dispersion structure for Toeplitz can (and often should) now be
    specified by \code{dispstrToep()}, see there.}
  \item{df}{integer value specifying the number of degrees of freedom
    of the multivariate t distribution used to construct the t copulas.}
  \item{df.fixed}{logical specifying if the degrees of freedom \code{df} will be
    considered as a parameter (to be estimated) or not.  The default,
    \code{FALSE}, means that \code{df} is to be estimated if the object is
    passed as argument to \code{\link{fitCopula}}.}
  \item{df.min}{non-negative number; the strict lower bound for
    \code{df}, mainly during fitting when \code{df.fixed=FALSE}, with
    \code{\link{fitCopula}}.}
  %% pCopula() methods:
  \item{copula}{an \R object of class \code{"\linkS4class{Copula}"}, in our
    case inheriting from \code{"\linkS4class{ellipCopula}"}.}

  \item{u}{a vector of the copula dimension \eqn{d} or a matrix with \eqn{d}
    columns, giving the points where the distribution function needs to be
    evaluated.  Note that values outside of the cube \eqn{[0,1]^d} are
    treated equivalently to those on the cube boundary.}
  \item{algorithm}{\code{NULL} or an \code{"algorithm"} object for package
    \CRANpkg{mvtnorm}'s \code{pmvt()} or \code{pmvnorm()} functions, see
    \code{\link[mvtnorm]{algorithms}}.  Note that for larger dimensions,
    the monte-carlo based \code{GenzBretz(..)} must be used, consequently
    with slightly random results.  By default, \code{algorithm = NULL},
    \code{algorithm} is chosen separately for each row \code{x <- u[i,]}, for
    \describe{
      \item{\code{normalCopula}:}{ via hidden function
	\code{pmvnormAlgo(dim, x, ...)}
	which currently is defined as % ../R/normalCopula.R
	\preformatted{ pmvnormAlgo <- function(dim, x, ...) {
    if(dim <= 3 && !anyNA(x) && (!any(xI <- x == Inf) || all(xI)))
        TVPACK(...)
    else if(dim <= 5)
        Miwa(...)
    else
        GenzBretz(...)
  }}
     }
      \item{\code{tCopula}:}{ via (hidden) \code{pmvtAlgo(dim, x, ...)}
	which is the same as \code{pmvnormAlgo()} above, but as \code{Miwa()}
	is not applicable, without the \code{else if(dim <= 5) Miwa(...)}
	clause.}
    }
  }
  \item{keepAttr}{\code{\link{logical}} passed to
    \code{\link[mvtnorm]{pmvnorm}} or \code{\link[mvtnorm]{pmvt}}, respectively.}
  \item{\dots}{for the \code{pCopula()} methods, optional arguments to the
    corresponding algorithm.}
  %% dispstrToep():
  \item{perm}{an \code{\link{integer}} vector of length \eqn{d =}\code{dim},
    which must be a permutation of \code{1:d} specifying the (column)
    ordering of the variables which has Toeplitz dispersion.}
  \item{check}{a \code{\link{logical}} specifying if the validity of
    \code{perm} should be checked (not strictly, currently).}
}
\value{For
  \describe{
    \item{\code{ellipCopula()}, \code{normalCopula()}, or \code{tCopula()}:}{
      an elliptical copula object of class \code{"\linkS4class{normalCopula}"}
      or \code{"\linkS4class{tCopula}"}.}
    \item{dispstrToep():}{the \code{\link{character}} string \code{"toep"},
      optionally with attribute (see \code{\link{attributes}},
      \code{\link{attr}}) \code{"perm"} with a permutation \eqn{p} of
      \code{1:d}, such that (the column permuted) \eqn{U_p}{U[p]}, or in
      the data case the column-permuted matrix \code{U[,p]} has as
      dispersion matrix \code{toeplitz(c(1, par))}, with \code{par} the
      respective parameter vector of bivariate \dQuote{correlations}
      \eqn{\rho_{i,j}}.

      Note that the result of \code{dispstrToep()} is currently stored in the
      \code{dispstr} slot of the copula object.}
    \item{pCopula(u, *): }{the numerical vector of copula probabilites of
      length \code{nrow(u)}.}
  }
}
\note{
  \code{ellipCopula()} is a wrapper for \code{normalCopula()} and
  \code{tCopula()}.

  The  \code{\link{pCopula}()} methods for the normal- and t-copulas
  accept optional arguments to be passed to the underlying
  (numerical integration) algorithms from package \CRANpkg{mvtnorm}'s
  \code{\link[mvtnorm]{pmvnorm}} and \code{\link[mvtnorm]{pmvt}},
  respectively, notably \code{algorithm}, see
  \code{\link[mvtnorm]{GenzBretz}}, or \code{abseps}
  which defaults to \code{0.001}.
}
%\author{Ivan Kojadinovic and Jun Yan}
\seealso{
  \code{\link{p2P}()}, and \code{\link{getSigma}()} for construction and
  extraction of the dispersion matrix \eqn{P} or \eqn{Sigma} matrix of
  (generalized)% for t
  correlations.

  \code{\link{archmCopula}}, \code{\link{fitCopula}}.
}
\examples{
normalCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un")
t.cop <- tCopula(c(0.5, 0.3), dim = 3, dispstr = "toep",
                 df = 2, df.fixed = TRUE)
getSigma(t.cop) # P matrix (with diagonal = 1)
stopifnot(all.equal(toeplitz(c(1, .5, .3)), getSigma(t.cop)))

## dispersion "AR1" :
nC.7 <- normalCopula(0.8, dim = 7, dispstr = "ar1")
getSigma(nC.7)
stopifnot(all.equal(toeplitz(.8^(0:6)), getSigma(nC.7)))

## from the wrapper
norm.cop <- ellipCopula("normal", param = c(0.5, 0.6, 0.7),
                        dim = 3, dispstr = "un")
if(require("scatterplot3d") && dev.interactive(orNone=TRUE)) {
  ## 3d scatter plot of 1000 random observations
  scatterplot3d(rCopula(1000, norm.cop))
  scatterplot3d(rCopula(1000, t.cop))
}
set.seed(12)
uN <- rCopula(512, norm.cop)
set.seed(2); pN1 <- pCopula(uN, norm.cop)
set.seed(3); pN2 <- pCopula(uN, norm.cop)
stopifnot(identical(pN1, pN2)) # no longer random for dim = 3
(Xtras <- copula:::doExtras())
if(Xtras) { ## a bit more accurately:
  set.seed(4); pN1. <- pCopula(uN, norm.cop, abseps = 1e-9)
  set.seed(5); pN2. <- pCopula(uN, norm.cop, abseps = 1e-9)
  stopifnot(all.equal(pN1., pN2., 1e-5))# see 3.397e-6
  ## but increasing the required precision (e.g., abseps=1e-15) does *NOT* help
}

## For smaller copula dimension 'd', alternatives are available and
## non-random, see ?GenzBretz from package 'mvtnorm' :
has_mvtn <- "package:mvtnorm" \%in\% search() #% << (workaround ESS Rd render bug)
if(!has_mvtn)
  require("mvtnorm")# -> GenzBretz(), Miva(), and TVPACK() are available
## Note that Miwa() would become very slow for dimensions 5, 6, ..
set.seed(4); pN1.M <- pCopula(uN, norm.cop, algorithm = Miwa(steps = 512))
set.seed(5); pN2.M <- pCopula(uN, norm.cop, algorithm = Miwa(steps = 512))
stopifnot(all.equal(pN1.M, pN2.M, tol= 1e-15))# *no* randomness
set.seed(4); pN1.T <- pCopula(uN, norm.cop, algorithm = TVPACK(abseps = 1e-10))
set.seed(5); pN2.T <- pCopula(uN, norm.cop, algorithm = TVPACK(abseps = 1e-14))
stopifnot(all.equal(pN1.T, pN2.T, tol= 1e-15))# *no* randomness (but no effect of 'abseps')
if(!has_mvtn)
   detach("package:mvtnorm")# (revert)


## Versions with unspecified parameters:
tCopula()
allEQ <- function(u,v) all.equal(u, v, tolerance=0)
stopifnot(allEQ(ellipCopula("norm"), normalCopula()),
          allEQ(ellipCopula("t"), tCopula()))
tCopula(dim=3)
tCopula(dim=4, df.fixed=TRUE)
tCopula(dim=5, disp = "toep", df.fixed=TRUE)
normalCopula(dim=4, disp = "un")

## Toeplitz after *permutation* dispersions (new in copula 1.1-0) ---------
tpar <- c(7,5,3)/8 # *gives* pos.def.:
(ev <- eigen(toeplitz(c(1, tpar)), symmetric=TRUE, only.values=TRUE)$values)
stopifnot(ev > 0)
N4.   <- ellipCopula("normal", dim=4, param=tpar, dispstr = "toep") #"regular"
## reversed order is "the same" for toeplitz structure:
N4.pr <- ellipCopula("normal", dim=4, param=tpar, dispstr = dispstrToep(4:1))
N4.p1 <- ellipCopula("normal", dim=4, param=tpar, dispstr = dispstrToep(c(4,1:3)))
N4.p2 <- ellipCopula("normal", dim=4, param=tpar, dispstr = dispstrToep(c(4:3,1:2)))
N4.p3 <- ellipCopula("normal", dim=4, param=tpar, dispstr = dispstrToep(c(2,4,1,3)))

(pm <- attr(N4.p3@dispstr, "perm")) # (2 4 1 3)
ip <- c(3,1,4,2) # the *inverse* permutation of (2 4 1 3) = Matrix::invPerm(pm)
(Sp3 <- getSigma(N4.p3)) # <-- "permuted toeplitz"
Sp3[ip, ip] # re-ordered rows & columns => *is* toeplitz :
stopifnot(exprs = {
  ## permutations  pm  and  ip  are inverses:
  pm[ip] == 1:4
  ip[pm] == 1:4
  is.matrix(T4 <- toeplitz(c(1, tpar)))
  identical(getSigma(N4.),   T4)
  identical(getSigma(N4.pr), T4) # 4:1 and 1:4 is "the same" for Rho
  identical(Sp3[ip, ip]  ,   T4)
  identical(Sp3,      T4[pm,pm])
})
## Data generation -- NB: The U matrices are equal only "in distribution":
set.seed(7); U.p3 <- rCopula(1000, N4.p3)
set.seed(7); U.   <- rCopula(1000, N4.)
stopifnot(exprs = {
 all.equal(loglikCopula(tpar, u=U.p3,      copula= N4.p3),
           loglikCopula(tpar, u=U.p3[,ip], copula= N4.) -> LL3)
 all.equal(loglikCopula(tpar, u=U.,      copula= N4.),
           loglikCopula(tpar, u=U.[,pm], copula= N4.p3) -> LL.)
})
c(LL. , LL3)# similar but different
if(Xtras) {
  fm.  <- fitCopula(N4.  , U.  )
  fm.3 <- fitCopula(N4.p3, U.p3)
  summary(fm.3)
  stopifnot(all.equal(coef(fm.), coef(fm.3), tol = 0.01))# similar but different
}
}
\keyword{distribution}
\keyword{multivariate}
