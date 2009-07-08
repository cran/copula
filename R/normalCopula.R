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


normalCopula <- function(param, dim = 2, dispstr = "ex") {
  pdim <- length(param)
  val <- new("normalCopula",
             dispstr = dispstr,
             dimension = dim,
             parameters = param,
             param.names = paste("rho", 1:pdim, sep="."),
             param.lowbnd = rep(-1, pdim),
             param.upbnd = rep(1, pdim),
             message = "Normal copula family",
             getRho = function(obj) {obj@parameters}
             )
  val
}


rnormalCopula <- function(copula, n) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  pnorm(rmvnorm(n, sigma = sigma))
}


pnormalCopula <- function(copula, u) {
  mycdf.vector <- function(x) {
    pmvnorm(lower = rep(-Inf, dim), upper = qnorm(x), sigma = sigma)
  }

  dim <- copula@dimension
  sigma <- getSigma(copula)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  u[u <= 0] <- 0
  u[u >= 1] <- 1
  val <- apply(u, 1, mycdf.vector)
  val
}

dnormalCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- qnorm(u)
  val <- dmvnorm(x, sigma = sigma) / apply(x, 1, function(v) prod(dnorm(v)))
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}


showNormalCopula <- function(object) {
  showCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
}


tailIndexNormalCopula <- function(copula) {
  rho <- copula@parameters
  upper <- lower <- ifelse(rho == 1, 1, 0)
  c(lower=lower, upper=upper)
}


kendallsTauNormalCopula <- function(copula) {
  rho <- copula@parameters
  2 * asin(rho) /pi
}

spearmansRhoNormalCopula <- function(copula) {
  rho <- copula@parameters
  asin(rho / 2) * 6 / pi
}

setMethod("rcopula", signature("normalCopula"), rnormalCopula)
setMethod("pcopula", signature("normalCopula"), pnormalCopula)
setMethod("dcopula", signature("normalCopula"), dnormalCopula)

setMethod("show", signature("normalCopula"), showNormalCopula)

setMethod("kendallsTau", signature("normalCopula"), kendallsTauNormalCopula)
setMethod("spearmansRho", signature("normalCopula"), spearmansRhoNormalCopula)
setMethod("tailIndex", signature("normalCopula"), tailIndexNormalCopula)

setMethod("calibKendallsTau", signature("normalCopula"), calibKendallsTauEllipCopula)
setMethod("calibSpearmansRho", signature("normalCopula"), calibSpearmansRhoEllipCopula)
