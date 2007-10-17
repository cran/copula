#################################################################
##   Copula R package by Jun Yan Copyright (C) 2007
##
##   This program is free software; you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation; either version 2 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program; if not, write to the Free Software Foundation, Inc.,
##   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
##
#################################################################

setClass("huslerReissCopula",
         representation = representation("evCopula"),
         contains = list("copula", "evCopula")
         )

AfunHuslerReiss <- function(copula, w) {
  alpha <- copula@parameters[1]
  w * pnorm(1 / alpha + 0.5 * alpha * log(w /(1 - w))) +
    (1 - w) * pnorm(1 / alpha - 0.5 * alpha * log(w / (1 - w))) 
}

huslerReissCopula <- function(param) {
  ## dim = 2
  dim <- 2
  val <- new("huslerReissCopula",
             dimension = dim,
             parameters = param[1],
             param.names = "param",
             param.lowbnd = 0,
             param.upbnd = Inf,
             message = "Husler-Reiss copula family; Extreme value copula")
  val
}


phuslerReissCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  ## Joe (1997, p.142)
  u1p <- -log(u1); u2p <- -log(u2)
  exp(- u1p * pnorm(1/alpha + 0.5 * log(u1p / u2p))
      - u2p * pnorm(1/alpha + 0.5 * log(u2p / u1p)))
}

dhuslerReissCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  ## Joe (1997, p.142)
  u1p <- -log(u1); u2p <- -log(u2); z <- u1p / u2p
  val <- 1/u1 / u2 * pcopula(copula, u) * 
    (pnorm(1/alpha - 0.5 * alpha * log(z)) *
     pnorm(1/alpha + 0.5 * alpha * log(z)) +
     0.5 * alpha / u2p * dnorm(1/alpha + 0.5 * alpha * log(z)))
  val
}


rhuslerReissCopula <- function(copula, n) {
  u1 <- runif(n)
  v <- runif(n)
  alpha <- copula@parameters[1]
  eps <- .Machine$double.eps ^ 0.8  ## don't know a better way
  myfun <- function(u2, u1, v) {
    ## Joe (1997, p.147)
    phuslerReissCopula(copula, cbind(u1, u2)) / u1 * pnorm(1/alpha + 0.5 * log(log(u1) / log(u2))) - v
  }
  ## I don't understand the rejection method used by finmetrics yet, so
  u2 <- sapply(1:n, function(x) uniroot(myfun, c(eps, 1 - eps), v=v[x], u1=u1[x])$root)
  cbind(u1, u2)
}


setMethod("pcopula", signature("huslerReissCopula"), phuslerReissCopula)
setMethod("dcopula", signature("huslerReissCopula"), dhuslerReissCopula)
setMethod("rcopula", signature("huslerReissCopula"), rhuslerReissCopula)

setMethod("Afun", signature("huslerReissCopula"), AfunHuslerReiss)
