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

tCopula <- function(param, dim = 2, dispstr = "ex", df = 4, df.fixed = FALSE) {
  pdim <- length(param)
  parameters <- param
  param.names <- paste("rho", 1:pdim, sep=".")
  param.lowbnd <- rep(-1, pdim)
  param.upbnd <- rep(1, pdim)
  if (!df.fixed) {
    parameters <- c(parameters, df)
    param.names <- c(param.names, "df")
    param.lowbnd <- c(param.lowbnd, 1) ## qt won't work for df < 1
    param.upbnd <- c(param.upbnd, Inf)
  }
  
  val <- new("tCopula",
             dispstr = dispstr,
             dimension = dim,
             parameters = parameters,
             df = df,
             df.fixed = df.fixed,
             param.names = param.names,
             param.lowbnd = param.lowbnd,
             param.upbnd = param.upbnd,
             message = paste("t copula family",
               if(df.fixed) paste("df fixed at", df) else NULL),
             getRho = function(obj) {
               if (df.fixed) obj@parameters
               else obj@parameters[-length(obj@parameters)]
             }
             )
  val
}

getdf <- function(object) {
  if (object@df.fixed) object@df
  else object@parameters[length(object@parameters)]
}

rtCopula <- function(copula, n) {
  dim <- copula@dimension
  df <- getdf(copula)
  sigma <- getSigma(copula)
  pt(rmvt(n, sigma = sigma, df = df), df = df)
}


ptCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- getdf(copula)
  mycdf.vector <- function(x) {
    pmvt(lower = rep(-Inf, dim), upper = qt(x, df = df), sigma = sigma, df = df)
  }
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  u[u <= 0] <- 0
  u[u >= 1] <- 1
  val <- apply(u, 1, mycdf.vector)
  val
}

dtCopula <- function(copula, u) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- getdf(copula)
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  x <- qt(u, df)
  val <- dmvt(x, delta = rep(0, dim), sigma = sigma, df = df, log = FALSE) /
    apply(x, 1, function(v) prod(dt(v, df = df)))
  val[apply(u, 1, function(v) any(v <= 0))] <- 0
  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}

showTCopula <- function(object) {
  showCopula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  if (object@df.fixed) cat("df is fixed at", object@df, "\n")
}

tailIndexTCopula <- function(copula) {
#### McNeil, Frey, Embrechts (2005), p.211
  df <- getdf(copula)
  rho <- copula@getRho(copula)
  upper <- lower <- 2 * pt(- sqrt((df + 1) * ( 1 - rho) / (1 + rho)), df=df + 1)
  c(upper=upper, lower=lower)
}

kendallsTauTCopula <- function(copula) {
  rho <- copula@getRho(copula)
  2 * asin(rho) /pi
}

spearmansRhoTCopula <- function(copula) {
  rho <- copula@getRho(copula)
  asin(rho / 2) * 6 / pi
}

setMethod("rcopula", signature("tCopula"), rtCopula)
setMethod("pcopula", signature("tCopula"), ptCopula)
setMethod("dcopula", signature("tCopula"), dtCopula)

setMethod("show", signature("tCopula"), showTCopula)

setMethod("kendallsTau", signature("tCopula"), kendallsTauTCopula)
setMethod("spearmansRho", signature("tCopula"), spearmansRhoTCopula)
setMethod("tailIndex", signature("tCopula"), tailIndexTCopula)

setMethod("calibKendallsTau", signature("tCopula"), calibKendallsTauEllipCopula)
setMethod("calibSpearmansRho", signature("tCopula"), calibSpearmansRhoEllipCopula)
