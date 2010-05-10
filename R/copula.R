#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

#### constructor
#### this function is outdated; may even be removed; ev in not covered yet
copula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("normal", "t", "clayton", "frank", "gumbel", "amh")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  if (fam <= 2) ellipCopula(family, param, dim, ...)
  else archmCopula(family, param, dim, ...)
}

#### show method
showCopula <- function(object) {
  cat(object@message, "\n")
  cat("Dimension: ", object@dimension, "\n")
  if (length(object@parameters) > 0) {
    cat("Parameters:\n")
    for (i in (1:length(object@parameters)))
      cat("  ", object@param.names[i], " = ", object@parameters[i], "\n")
  }
}

setMethod("show", signature("copula"), showCopula)


#### numerical computation of association measures
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


#### numerical tail index, not accurate
tailIndexCopula <- function(copula, eps = .Machine$double.eps^0.5) {
  u <- eps
  v <- 1 - u
  lower <- pcopula(copula, c(u, u))/u
  upper <- (1 - 2 * v + pcopula(copula, c(v, v)))/ (1 - v)
  c(lower=lower, upper=upper)
}

# setMethod("kendallsTau", signature("copula"), kendallsTauCopula)
# setMethod("spearmansRho", signature("copula"), spearmansRhoCopula)
setMethod("tailIndex", signature("copula"), tailIndexCopula)

#### numerical calibration
calibKendallsTauCopula <- function(copula, tau) {
  myfun <- function(theta) {
    copula@parameters <- theta
    kendallsTau(copula) - tau
  }
  .eps <- .Machine$double.eps^.5
  lower <- copula@param.lowbnd + .eps
  upper <- copula@param.upbnd - .eps
##   lower <- ifelse(lower == -Inf, -.Machine$double.xmax^.5, lower)
##   upper <- ifelse(upper == Inf, .Machine$double.xmax^.5, upper)
  lower <- ifelse(lower == -Inf, -sqrt(.Machine$double.xmax), lower)
  upper <- ifelse(upper == Inf, sqrt(.Machine$double.xmax), upper)
  sol <- uniroot(myfun, interval=c(lower, upper))$root
  sol
}

calibSpearmansRhoCopula <- function(copula, rho) {
  myfun <- function(theta) {
    copula@parameters <- theta
    spearmansRho(copula) - rho
  }
  .eps <- .Machine$double.eps^.5
  lower <- copula@param.lowbnd
  upper <- copula@param.upbnd
##   lower <- ifelse(lower == -Inf, -.Machine$double.xmax^.5, lower)
##   upper <- ifelse(upper == Inf, .Machine$double.xmax^.5, upper)
  lower <- ifelse(lower == -Inf, -sqrt(.Machine$double.xmax), lower)
  upper <- ifelse(upper == Inf, sqrt(.Machine$double.xmax), upper)
  sol <- uniroot(myfun, interval=c(lower, upper))$root
  sol
}

# setMethod("calibKendallsTau", signature("copula"), calibKendallsTauCopula)
# setMethod("calibSpearmansRho", signature("copula"), calibSpearmansRhoCopula)
