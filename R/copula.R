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


##' show method
showCopula <- function(object) {
  cat(object@message, "\n")
  cat("Dimension: ", object@dimension, "\n")
  if (length(object@parameters) > 0) {
    cat("Parameters:\n")
    for (i in (1:length(object@parameters)))
      cat("  ", object@param.names[i], " = ", object@parameters[i], "\n")
  }
  invisible(object)
}

setMethod("show", signature("copula"), showCopula)


### numerical computation of association measures

## kendallsTauCopula <- function(copula, eps = NULL, ...) {
##   integrand <- function(u) pcopula(copula, u) * dcopula(copula, u)
##   if (is.null(eps)) .eps <- .Machine$double.eps^0.9
##   else .eps <- eps
##   lower <- c(.eps, .eps)
##   upper <- c(1 - .eps, 1 - .eps)
##   integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
##   4 * integ - 1
## }

## spearmansRhoCopula <- function(copula, eps = NULL, ...) {
##   integrand <- function(u) pcopula(copula, u)
##   if (is.null(eps)) .eps <- .Machine$double.eps^0.9
##   else .eps <- eps
##   lower <- c(.eps, .eps)
##   upper <- c(1 - .eps, 1 - .eps)
##   integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
##   12 * integ - 3
## }


### numerical tail index, not accurate

tailIndexCopula <- function(copula, eps = .Machine$double.eps^0.5) {
  u <- eps
  v <- 1 - u
  lower <- pcopula(copula, c(u, u))/u
  upper <- (1 - 2 * v + pcopula(copula, c(v, v)))/ u
  c(lower=lower, upper=upper)
}

# setMethod("kendallsTau", signature("copula"), kendallsTauCopula)
# setMethod("spearmansRho", signature("copula"), spearmansRhoCopula)
setMethod("tailIndex", signature("copula"), tailIndexCopula)


### numerical calibration

calibKendallsTauCopula <- function(copula, tau) {
  myfun <- function(theta) {
    copula@parameters <- theta
    kendallsTau(copula) - tau
  }
  .eps <- .Machine$double.eps^.5
  lower <- pmax(-sqrt(.Machine$double.xmax), copula@param.lowbnd + .eps)
  upper <- pmin( sqrt(.Machine$double.xmax), copula@param.upbnd  - .eps)
  uniroot(myfun, interval=c(lower, upper))$root
}

calibSpearmansRhoCopula <- function(copula, rho) {
  myfun <- function(theta) {
    copula@parameters <- theta
    spearmansRho(copula) - rho
  }
  lower <- pmax(-sqrt(.Machine$double.xmax), copula@param.lowbnd)
  upper <- pmin( sqrt(.Machine$double.xmax), copula@param.upbnd )
  uniroot(myfun, interval=c(lower, upper))$root
}

# setMethod("calibKendallsTau", signature("copula"), calibKendallsTauCopula)
# setMethod("calibSpearmansRho", signature("copula"), calibSpearmansRhoCopula)

###-- "Copula" methods + glue  former "copula" <--> former "nacopula" ---------

setMethod("dim", "copula",
	  function(x) x@dimension)

## Dummy bail-out methods for all generics --> ./zzz.R
##  "nacopula" methods                     --> ./nacopula.R
