## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


##' print() / show() method, _or_ ("tCopula", e.g.) an auxiliary for that :
printCopula <- function(x, digits = getOption("digits"),
                        descr.kind = c("short", "very short", "long"), ...) {
  validObject(x)
  descr.kind <- match.arg(descr.kind)
  cat(describeCop(x, kind = descr.kind), "\n")
  cat("Dimension: ", dim(x), "\n")
  if (length(par <- getTheta(x, freeOnly=FALSE, attr = TRUE)) > 0) {
    ## hasFx <- !is.null(.fixed <- attr(par, "fixed")) && any(.fixed)
    ## JY: no attribute "fixed" is available from getTheta
    .fixed <- !isFree(x); hasFx <- any(.fixed)
    hasFx <- any(!isFree(x))
    cat(sprintf("Parameters%s:\n",
		if(hasFx) " (partly fixed, see ':=')" else ""))
    ## FIXME:  for the d=400  ellipsCopula with dispstr = "un": do *not* print all!
    pnms <- format(names(par)) # padding
    pars <- format(as.vector(par), digits=digits)
    for (i in seq_along(par))
      cat(sprintf("  %s %3s %s\n",
		  pnms[i], if(hasFx && .fixed[i]) ":=" else " =",
		  pars[i]))
  }
  invisible(x)
}

## default print() and show() method for copulas:
print.parCopula <- printCopula
setMethod("show", "parCopula", function(object) printCopula(object))


### numerical computation of association measures

## tauCopula <- function(copula, eps = NULL, ...) {
##   integrand <- function(u) pCopula(u, copula) * dCopula(u, copula)
##   if (is.null(eps)) .eps <- .Machine$double.eps^0.9
##   else .eps <- eps
##   lower <- c(.eps, .eps)
##   upper <- c(1 - .eps, 1 - .eps)
##   integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
##   4 * integ - 1
## }

## rhoCopula <- function(copula, eps = NULL, ...) {
##   integrand <- function(u) pCopula(u, copula)
##   if (is.null(eps)) .eps <- .Machine$double.eps^0.9
##   else .eps <- eps
##   lower <- c(.eps, .eps)
##   upper <- c(1 - .eps, 1 - .eps)
##   integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
##   12 * integ - 3
## }


### numerical tail index, not accurate
## lambdaCopula <- function(copula, eps = .Machine$double.eps^0.5) {
##     u <- eps
##     v <- 1 - u
##     lower <- pCopula(c(u, u), copula)/u
##     upper <- (1 - 2 * v + pCopula(c(v, v), copula))/ u
##     c(lower=lower, upper=upper)
## }

## setMethod("tau", signature("copula"), tauCopula)
## setMethod("rho", signature("copula"), rhoCopula)
## setMethod("lambda", signature("copula"), lambdaCopula)


### Numerical "calibration": inverse tau() and rho()

## uniroot()'s  precision tol = .Machine$double.eps^.25 = 0.000122 ~= 1e-4  is a bit too small (for MM)
## -> decreasing to   tol = 1e-7 [2014-05-20]
iTauCopula <- function(copula, tau, bound.eps = .Machine$double.eps^.5, tol = 1e-7, ...) {
  stopifnot(length(bound.eps) == 1, is.finite(bound.eps), 0 <= bound.eps,
	    2*bound.eps < copula@param.upbnd - copula@param.lowbnd)
  myfun <- function(theta) {
    copula@parameters <- theta
    tau(copula) - tau
  }
  lower <- pmax(-sqrt(.Machine$double.xmax), copula@param.lowbnd + bound.eps)
  upper <- pmin( sqrt(.Machine$double.xmax), copula@param.upbnd  - bound.eps)
  uniroot(myfun, interval=c(lower, upper), tol=tol, ...)$root
}

iRhoCopula <- function(copula, rho, bound.eps = 0, tol = 1e-7, ...) {
  stopifnot(length(bound.eps) == 1, is.finite(bound.eps), 0 <= bound.eps,
	    2*bound.eps < copula@param.upbnd - copula@param.lowbnd)
  d_rho <- function(theta) {
    copula@parameters <- theta
    rho(copula) - rho
  }
  lower <- pmax(-sqrt(.Machine$double.xmax), copula@param.lowbnd + bound.eps)
  upper <- pmin( sqrt(.Machine$double.xmax), copula@param.upbnd  - bound.eps)
  uniroot(d_rho, interval = c(lower, upper), tol=tol, ...)$root
}

## From Ivan: do not make these the default methods for *all* copulas, as
## some have more than one parameter
setMethod("iTau", signature("archmCopula"), iTauCopula)
setMethod("iRho", signature("archmCopula"), iRhoCopula)


###-- "Copula" methods + glue  former "copula" <--> former "nacopula" ---------

## Dummy bail-out methods for all generics --> ./zzz.R
##  "nacopula" methods                     --> ./nacopula.R

setGeneric("dPsi", function(copula, ...) standardGeneric("dPsi"))
setMethod("dPsi", "acopula",
          function(copula, t, theta, degree=1, log=FALSE, ...) {
	      s <- if(log || degree %% 2 == 0) 1. else -1.
	      s * copula@absdPsi(t, theta, degree=degree, log=log, ...)
       })

## Methods for 'xcopula': extended copulas defined by mainly *one* copula slot
setMethod("dim", signature("xcopula"), function(x) dim(x@copula))

## logical indicating which parameters are free
setMethod("isFree", signature("rotCopula"), function(copula) isFree(copula@copula))

## number of (free / all) parameters :
setMethod("nParam", signature("rotCopula"), function(copula, freeOnly=FALSE)
    nParam(copula@copula, freeOnly=freeOnly))

## paramNames() methods: further provided separately for "fitCopula", "fitMvdc"
setMethod("paramNames", signature("copula"), function(x) x@param.names[isFreeP(x@parameters)])


