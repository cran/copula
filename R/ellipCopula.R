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

getSigma <- function(copula) {
  dim <- copula@dimension
  rho <- copula@getRho(copula)
  sigma <- diag(dim)
  if (copula@dispstr == "ex") {
    sigma[lower.tri(sigma)] <- rho[1]
    sigma[upper.tri(sigma)] <- rho[1]
  }
  else if (copula@dispstr == "ar1") {
    for (i in 1:dim)  for (j in 1:dim)  sigma[i,j] <- rho ^ abs(i - j)
  }
  else if (copula@dispstr == "un") {
    sigma[lower.tri(sigma)] <- rho
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  }
  else if (copula@dispstr == "toep") {
    for (i in 1:dim) for (j in 1:dim)
      if (i != j) sigma[i,j] <- rho[abs(i - j)]
  }
  sigma
}

ellipCopula <- function(family, param, dim = 2, dispstr = "ex", df = 4, ...) {
  familiesImplemented <- c("normal", "t")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   normalCopula(param, dim = dim, dispstr = dispstr),
                   tCopula(param, dim = dim, dispstr = dispstr, df = df, ...)
                   )
  copula
}

calibKendallsTauEllipCopula <- function(copula, tau) {
  sin((tau * pi) / 2)
}

calibSpearmansRhoEllipCopula <- function(copula, rho) {
  sin (pi * rho / 6) * 2
}

tauDerEllipCopula <- function(copula)  {
  return( 2 / (pi * sqrt(1 - copula@parameters^2)) )
}

tauDerFunEllipCopula <- function(copula)  {
  function(x) return( 2 / (pi * sqrt(1 - x^2)) )
}

rhoDerEllipCopula <- function(copula) {
  return( 6 / (pi * sqrt(4 - copula@parameters^2)) )
}

rhoDerFunEllipCopula <- function(copula) {
  function(x) return( 6 / (pi * sqrt(4 - x^2)) )
}

setMethod("calibKendallsTau", signature("ellipCopula"), calibKendallsTauEllipCopula)
setMethod("calibSpearmansRho", signature("ellipCopula"), calibSpearmansRhoEllipCopula)

setMethod("tauDer", signature("ellipCopula"), tauDerEllipCopula)
setMethod("rhoDer", signature("ellipCopula"), rhoDerEllipCopula)

setMethod("tauDerFun", signature("ellipCopula"), tauDerFunEllipCopula)
setMethod("rhoDerFun", signature("ellipCopula"), rhoDerFunEllipCopula)
