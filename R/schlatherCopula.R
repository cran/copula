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


## Schlather copula; does not offer full range of dependence
setClass("schlatherCopula",
         representation = representation("evCopula"),
           # exprdist = "expression"),
         contains = list("evCopula")
         )

AfunSchlather <- function(copula, w) { ## one-parameter for now
  alpha <- copula@parameters[1]
  A <- 0.5 * (1 + sqrt(1 - 2 * (alpha + 1) * w * (1 - w)))
  ifelse(w == 0 | w == 1, 1, A)
}

AfunDerSchlather <- function(copula, w) {
  alpha <- copula@parameters[1]
  ainv <- 1 / alpha
  z <- 0.5 * alpha * log(w / (1 - w))
  ## deriv(~ 0.5 * (1 + (1 - 2 * (alpha + 1) * w * (1 - w))), "w", hessian=TRUE)
  zder <- eval(expression({
    .expr2 <- 2 * (alpha + 1)
    .expr3 <- .expr2 * w
    .expr4 <- 1 - w
    .value <- 0.5 * (1 + (1 - .expr3 * .expr4))
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("w"), c("w")))
    .grad[, "w"] <- -(0.5 * (.expr2 * .expr4 - .expr3))
    .hessian[, "w", "w"] <- 0.5 * (.expr2 + .expr2)
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha = alpha, w = w))
  der1 <- c(attr(zder, "gradient"))
  der2 <- c(attr(zder, "hessian"))
  data.frame(der1 = der1, der2 = der2)
}

derAfunWrtParamSchlather <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## to be completed
  stop("to be implemented")
}

schlatherCopula <- function(param) {
  dim <- 2L
  new("huslerReissCopula",
             dimension = dim,
             ## exprdist = c(cdf = cdf, pdf = pdf),
             parameters = param[1],
             param.names = "param",
             param.lowbnd = -1,
             param.upbnd = 1,
             message = "Schlather copula family; Extreme value copula")
}


pschlatherCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  ## Beirlant, Goegebeur, Segers, and Teugels (2004, p.295)
  w <- log(u2) / log(u1 * u2)
  u1 * u2 * exp(AfunSchlather(copula, w))
}

dschlatherCopula <- function(copula, u, log=FALSE, ...) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  alpha <- copula@parameters[1]
  stop("to be implemented")
}


## This block is copied from ../../copulaUtils/assoc/ ##########################

#setMethod("pcopula", signature("schlatherCopula"), pschlatherCopula)
#setMethod("dcopula", signature("schlatherCopula"), dschlatherCopula)
## revCopula is much faster
## setMethod("rcopula", signature("schlatherCopula"), rschlatherCopula)

#setMethod("Afun", signature("schlatherCopula"), AfunSchlather)
#setMethod("AfunDer", signature("schlatherCopula"), AfunDerSchlather)

## setMethod("kendallsTau", signature("schlatherCopula"), kendallsTauHuslerReissCopula)
## setMethod("spearmansRho", signature("schlatherCopula"), spearmansRhoHuslerReissCopula)

## setMethod("calibKendallsTau", signature("schlatherCopula"), calibKendallsTauHuslerReissCopula)
## setMethod("calibSpearmansRho", signature("schlatherCopula"), calibSpearmansRhoHuslerReissCopula)

## setMethod("tauDer", signature("schlatherCopula"), tauDerHuslerReissCopula)
## setMethod("rhoDer", signature("schlatherCopula"), rhoDerHuslerReissCopula)
