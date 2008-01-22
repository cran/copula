#################################################################################
##
##   R package Copula by Jun Yan Copyright (C) 2008
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


perspCopula <- function(x, fun, n = 51, theta = -30, phi = 30, expand = 0.618, ...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}



contourCopula <- function(x, fun, n = 51,...) {
  eps <- (.Machine$double.eps)^(1/4)
  eps <- 0
  xis <- yis <- seq(0 + eps, 1 - eps, len = n)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), n, n)
  contour(xis, yis, zmat, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}


perspMvdc <- function(x, fun,
                      xlim, ylim, nx = 51, ny = 51,
                      theta = -30, phi = 30, expand = 0.618, ...) {
  xis <- seq(xlim[1], xlim[2], length = nx)
  yis <- seq(ylim[1], ylim[2], length = ny)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), nx, ny)
  persp(xis, yis, zmat, theta = theta, phi = phi, expand = expand, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}



contourMvdc <- function(x, fun,
                        xlim, ylim, nx = 51, ny = 51, ...) {
  xis <- seq(xlim[1], xlim[2], length = nx)
  yis <- seq(ylim[1], ylim[2], length = ny)
  grids <- as.matrix(expand.grid(xis, yis))
  zmat <- matrix(fun(x, grids), nx, ny)
  contour(xis, yis, zmat, ...)
  val <- list(x = xis, y = yis, z = zmat)
  invisible(val)
}

setMethod("persp", signature("copula"), perspCopula)
setMethod("contour", signature("copula"), contourCopula)

setMethod("persp", signature("mvdc"), perspMvdc)
setMethod("contour", signature("mvdc"), contourMvdc)
