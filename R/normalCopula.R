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


normalCopula <- function(param = NA_real_, dim = 2L, dispstr = "ex") {
  stopifnot((pdim <- length(param)) >= 1)
  new("normalCopula",
      dispstr = dispstr,
      dimension = as.integer(dim),
      parameters = param,
      param.names = paste("rho", 1:pdim, sep="."),
      param.lowbnd = rep(-1, pdim),
      param.upbnd = rep(1, pdim),
      fullname = "Normal copula family",
      getRho = function(obj) obj@parameters)
}


rnormalCopula <- function(n, copula)
    pnorm(rmvnorm(n, sigma = getSigma(copula)))


pnormalCopula <- function(u, copula) {
  dim <- copula@dimension
  i.lower <- rep.int(-Inf, dim)
  sigma <- getSigma(copula)
  u <- matrix(pmax(0, pmin(1, u)), ncol = dim)
  apply(qnorm(u), 1, function(qu)
        pmvnorm(lower = i.lower, upper = qu, sigma = sigma))
}

dnormalCopula <- function(u, copula, log=FALSE, ...) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  x <- qnorm(u)
  ## work in log-scale [less over-/under-flow, then (maybe) transform:
  val <- dmvnorm(x, sigma = sigma, log=TRUE) - rowSums(dnorm(x, log=TRUE))
  if(any(out <- !is.na(u) & (u <= 0 | u >= 1)))
    val[apply(out, 1, any)] <- -Inf
  if(log) val else exp(val)
}


showNormalCopula <- function(object) {
  showCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  invisible(object)
}


tailIndexNormalCopula <- function(copula) {
  rho <- copula@parameters
  upper <- lower <- ifelse(rho == 1, 1, 0)
  c(lower=lower, upper=upper)
}


tauNormalCopula <- function(copula) {
  rho <- copula@parameters
  2 * asin(rho) /pi
}

rhoNormalCopula <- function(copula) {
  rho <- copula@parameters
  asin(rho / 2) * 6 / pi
}

setMethod("rCopula", signature("numeric", "normalCopula"), rnormalCopula)

setMethod("pCopula", signature("matrix", "normalCopula"), pnormalCopula)
setMethod("pCopula", signature("numeric", "normalCopula"),pnormalCopula)
setMethod("dCopula", signature("matrix", "normalCopula"), dnormalCopula)
setMethod("dCopula", signature("numeric", "normalCopula"),dnormalCopula)

setMethod("show", signature("normalCopula"), showNormalCopula)

setMethod("tau", signature("normalCopula"), tauNormalCopula)
setMethod("rho", signature("normalCopula"), rhoNormalCopula)
setMethod("tailIndex", signature("normalCopula"), tailIndexNormalCopula)

setMethod("iTau", signature("normalCopula"), iTauEllipCopula)
setMethod("iRho", signature("normalCopula"), iRhoEllipCopula)
