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


AfunGalambos <- function(copula, w) {
  alpha <- copula@parameters[1]
  1 - (w^(-alpha) + (1 - w)^(-alpha))^(-1/alpha)
}

galambosCopula <- function(param) {
  ## dim = 2
  dim <- 2
  val <- new("galambosCopula",
             dimension = dim,
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             message = "Galambos copula family; Extreme value copula")
  val
}
  
pgalambosCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(galambosCopula.algr$cdf)
}

dgalambosCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(galambosCopula.algr$pdf)
}


rgalambosCopula <- function(copula, n) {
  u1 <- runif(n)
  v <- runif(n)
  alpha <- copula@parameters[1]
  deriv1 <- function(u1, u2) {
    eval(galambosCopula.algr$deriv1cdf)
  }
  eps <- .Machine$double.eps ^ 0.8  ## don't know a better way
  myfun <- function(u2, u1, v) {
    deriv1(u1, u2)/deriv1(u1, 1 - eps) - v
  }
  ## I don't understand the rejection method used by finmetrics yet, so
  u2 <- sapply(1:n, function(x) uniroot(myfun, c(eps, 1 - eps), v=v[x], u1=u1[x])$root)
  cbind(u1, u2)
}


setMethod("pcopula", signature("galambosCopula"), pgalambosCopula)
setMethod("dcopula", signature("galambosCopula"), dgalambosCopula)
setMethod("rcopula", signature("galambosCopula"), rgalambosCopula)

setMethod("Afun", signature("galambosCopula"), AfunGalambos)
