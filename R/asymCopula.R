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


### asymmetric copulas #########################################################

setClass("asymCopula", contains = "copula",
         representation = representation(
           copula1 = "copula",
           copula2 = "copula"
           ),
         validity = function(object) {
             if(object@copula1@dimension != object@copula2@dimension)
                 return("The dimensions of the two copulas are different")
           dim <- object@dimension
           ## from Classes.R
           if (dim != as.integer(dim))
             return("dim must be integer")
           if (dim < 2)
             return("dim must be >= 2")
           param <- object@parameters
           upper <- object@param.upbnd
           lower <- object@param.lowbnd
           if (length(param) != length(upper))
             return("Parameter and upper bound have non-equal length")
           if (length(param) != length(lower))
             return("Parameter and lower bound have non-equal length")
           if (any(is.na(param) | param > upper | param < lower))
             return("Parameter value out of bound")
	   ## else return
	   TRUE
	 })


## Liebscher (2008, JMA); the special case is Khoudraji
gfun <- function(u, a) u^a
igfun <- function(u, a) u^(1 / a)
gfunDer <- function(u, a) a * u^(a - 1)

asymCopula <- function(shapes, copula1, copula2) {
  val <- new("asymCopula",
             dimension = copula1@dimension,
             parameters = c(copula1@parameters, copula2@parameters, shapes),
             param.names = c(copula1@param.names, copula2@param.names, "shape1", "shape2"),
             param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, 0, 0),
             param.upbnd = c(copula1@param.upbnd, copula2@param.upbnd, 1, 1),
             copula1 = copula1,
             copula2 = copula2,
             message = "Asymmetric Copula")
}


getCopulaComps <- function(object) {
  p1 <- length(object@copula1@parameters)
  p2 <- length(object@copula2@parameters)
  shape <- object@parameters[(p1 + p2) + 1:2]
  copula1 <- object@copula1
  copula2 <- object@copula2
  if (p1 > 0) slot(copula1, "parameters") <- object@parameters[1:p1]
  if (p2 > 0) slot(copula2, "parameters") <- object@parameters[p1 + 1:p2]
  list(shape = shape, copula1 = copula1, copula2 = copula2)
}

AfunAsymCopula <- function(copula, w) {
  comps <- getCopulaComps(copula)
  copula1 <- comps$copula1; copula2 <- comps$copula2
  ## assuming copula1 and copula2 are both evCopula
  stopifnot(is(copula1, "evCopula"), is(copula2, "evCopula"))
  a1 <- comps$shape[1];  a2 <- comps$shape[2]
  den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
  den2 <- a1 * (1 - w) + a2 * w
  t1 <- (1 - a2) * w / den1; t1 <- ifelse(is.na(t1), 1, t1)
  t2 <- a2 * w / den2; t2 <- ifelse(is.na(t2), 1, t2)
  den1 * Afun(copula1, t1) + den2 * Afun(copula2, t2)
}

pasymCopula <- function(copula, u) {
  if(!is.matrix(u)) u <- matrix(u, ncol = 2)
  comps <- getCopulaComps(copula)
  a1 <- comps$shape[1];  a2 <- comps$shape[2]
  copula1 <- comps$copula1; copula2 <- comps$copula2
  gu1 <- cbind(gfun(u[,1], 1 - a1), gfun(u[,2], 1 - a2))
  gu2 <- cbind(gfun(u[,1], a1), gfun(u[,2], a2))
  pcopula(copula1, gu1) * pcopula(copula2, gu2)
}

dasymCopula <- function(copula, u, log=FALSE, ...) {
  ## WARNING:
  ## The following derivation assumes that both components are symmetric!
  ## See dC1du and dC2du; they don't distinguish u1 or u2.
  if(!is.matrix(u)) u <- matrix(u, ncol = 2)
  comps <- getCopulaComps(copula)
  if(log) stop("'log=TRUE' not yet implemented")
  a1 <- comps$shape[1];  a2 <- comps$shape[2]
  copula1 <- comps$copula1; copula2 <- comps$copula2
  gu1 <- cbind(gfun(u[,1], 1 - a1), gfun(u[,2], 1 - a2))
  gu2 <- cbind(gfun(u[,1], a1), gfun(u[,2], a2))
  dC1du <- derCdfWrtArgs(copula1, gu1)
  dC2du <- derCdfWrtArgs(copula2, gu2)
  part1 <- dcopula(copula1, gu1) * gfunDer(u[,1], 1 - a1) * gfunDer(u[,2], 1 - a2) * pcopula(copula2, gu2)
  part2 <- dC1du[,1] * gfunDer(u[,1], 1 - a1) * gfunDer(u[,2], a2) * dC2du[,2]
  part3 <- dC1du[,2] * gfunDer(u[,2], 1 - a2) * gfunDer(u[,1], a1) * dC2du[,1]
  part4 <- pcopula(copula1, gu1) * dcopula(copula2, gu2) * gfunDer(u[,2], a2) * gfunDer(u[,1], a1)
  part1 + part2 + part3 + part4
}

rasymCopula <- function(copula, n) {
  comps <- getCopulaComps(copula)
  a1 <- comps$shape[1];  a2 <- comps$shape[2]
  copula1 <- comps$copula1; copula2 <- comps$copula2
  ## Theorem 2.1, Lemma 2.1, Liebscher (2008, JMA)
  u <- rcopula(copula1, n)
  v <- rcopula(copula2, n)
  x <- matrix(NA, n, 2)
  x[,1] <- pmax(igfun(u[,1], 1 - a1), igfun(v[,1], a1))
  x[,2] <- pmax(igfun(u[,2], 1 - a2), igfun(v[,2], a2))
  x
}

setMethod("Afun", signature("asymCopula"), AfunAsymCopula)

setMethod("rcopula", signature("asymCopula"), rasymCopula)
setMethod("pcopula", signature("asymCopula"), pasymCopula)
setMethod("dcopula", signature("asymCopula"), dasymCopula)

