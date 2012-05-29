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


### asymmetric explicit copulas ################################################

setClass("asymExplicitCopula",
         representation = representation(
           "copula",
           copula1 = "copula",
           copula2 = "copula",
           exprdist = "expression",
           derExprs1 = "expression",
           derExprs2 = "expression"
           ),
         validity = function(object) {
           stopifnot(object@copula1@dimension == object@copula2@dimension)
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
           else return (TRUE)
         },
         contains = list("copula")
         )

## Liebscher (2008, JMA); the special case is Khoudraji
## gfun <- function(u, a) u^a
## igfun <- function(u, a) u^(1 / a)
## gfunDer <- function(u, a) a * u^(a - 1)

## only for power function in Liebscher (2008)

asymExplicitCopula <- function(shapes, copula1, copula2) {
  d <- copula1@dimension
  stopifnot(d == copula2@dimension,
            d == length(shapes))
  ## cdf
  getcdfchar <- function(cdf, om=FALSE) {
    ## WARNING: this only works up to dim 9; e.g., u10 could be replaced with u1^shp1
    if (d >= 10) stop("The maximum implemented dim is 9.")
    cdf <- deparse(cdf)
    for (i in 1:d) {
      ui <- paste("u", i, sep="")
      shpi <- paste("shp", i, sep="")
      if (om) shpi <- paste("(1 - ", shpi, ")")
      replacement <- paste("(", ui, "^", shpi, ")")
      cdf <- gsub(ui, replacement, cdf)
    }
    cdf
  }

  cdf1 <- getcdfchar(copula1@exprdist$cdf, TRUE)
  cdf2 <- getcdfchar(copula2@exprdist$cdf, FALSE)
  cdf <- c("(", cdf1, ") * (", cdf2, ")")
  cdf <- parse(text = cdf)
  ## pdf
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }
  if (d <= 6) pdf <- pdfExpr(cdf, d)
  else {
    warnings("The pdf is only available for dim 6 or lower.")
    pdf <- NULL
  }
  derExprs <- function(cdf, n) {
    val <- as.expression(cdf)
    for (i in 1:n) {
      val <- c(val, D(val[i], paste("u", i , sep="")))
    }
    val
  }
  derExprs1 <- derExprs(copula1@exprdist$cdf, d)
  derExprs2 <- derExprs(copula2@exprdist$cdf, d)
  shapes.names <- paste("shape", 1:d, sep="")
  val <- new("asymExplicitCopula",
             dimension = d,
             parameters = c(copula1@parameters, copula2@parameters, shapes),
             param.names = c(copula1@param.names, copula2@param.names, shapes.names),
             param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, rep(0, d)),
             param.upbnd = c(copula1@param.upbnd, copula2@param.upbnd, rep(1, d)),
             copula1 = copula1,
             copula2 = copula2,
             exprdist = c(cdf=cdf, pdf=pdf),
             derExprs1 = derExprs1, derExprs2 = derExprs2,
             message = "Asymmetric Explicit Copula")
}


getAsymExplicitCopulaComps <- function(object) {
  p1 <- length(object@copula1@parameters)
  p2 <- length(object@copula2@parameters)
  d <- object@dimension
  shapes <- object@parameters[(p1 + p2) + 1:d]
  copula1 <- object@copula1
  copula2 <- object@copula2
  if (p1 > 0) slot(copula1, "parameters") <- object@parameters[1:p1]
  if (p2 > 0) slot(copula2, "parameters") <- object@parameters[p1 + 1:p2]
  list(shapes = shapes, copula1 = copula1, copula2 = copula2)
}

## AfunAsymCopula <- function(copula, w) {
##   ## assuming copula@copula1 and copula@copula2 are both evCopula
##   comps <- getCopulaComps(copula)
##   a1 <- comps$shape[1];  a2 <- comps$shape[2]
##   copula1 <- comps$copula1; copula2 <- comps$copula2
##   den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
##   den2 <- a1 * (1 - w) + a2 * w
##   t1 <- (1 - a2) * w / den1; t1 <- ifelse(is.na(t1), 1, t1)
##   t2 <- a2 * w / den2; t2 <- ifelse(is.na(t2), 1, t2)
##   den1 * Afun(copula1, t1) + den2 * Afun(copula2, t2)
## }

pasymExplicitCopula <- function(copula, u) {
  u <- as.matrix(u)
  comps <- getAsymExplicitCopulaComps(copula)
  p1 <- pcopula(comps$copula1, t(t(u)^(1 - comps$shapes)))
  p2 <- pcopula(comps$copula2, t(t(u)^comps$shapes))
  p1 * p2
}

getPowerSet <- function(d) {
  TF <- matrix(c(TRUE, FALSE), 2, d)
  as.matrix(expand.grid(as.list(as.data.frame(TF))))
}

densDers <- function(idx, u, dg, copula, derExprs) {
  ## assuming exchangeable copula1 and copula2
  dorder <- sum(idx)
  alpha <- copula@parameters[1]
  d <- copula@dimension
  newidx <- c((1:d)[idx], (1:d)[!idx])
  u <- u[, newidx]
  for (i in 1:d) assign(paste("u", i, sep=""), u[,i])
  dgu <- if(sum(idx) == 0) 1 else apply(dg[,idx,drop=FALSE], 1, prod)
  c(eval(derExprs[dorder + 1])) * dgu
}

dasymExplicitCopula <- function(copula, u, log=FALSE, ...) {
  u <- as.matrix(u)
  comps <- getAsymExplicitCopulaComps(copula)
  a <- comps$shapes
  if(log) stop("'log=TRUE' not yet implemented")
  u1 <- t(t(u)^(1 - a))
  u2 <- t(t(u)^a)
  dg1 <- t((1 - a) * t(u)^(-a))
  dg2 <- t(a * t(u)^(a - 1))
  d <- copula@dimension
  powerSet <- getPowerSet(d)
  dens <- 0
  for (i in 1:nrow(powerSet)) {
    idx1 <- c(powerSet[i,])
    part1 <- densDers(idx1, u1, dg1, copula@copula1, copula@derExprs1)
    idx2 <- c(!powerSet[i,])
    part2 <- densDers(idx2, u2, dg2, copula@copula2, copula@derExprs2)
    dens <- dens + part1 * part2
    ## print(part1); print(part2)
  }
  dens
}

rasymExplicitCopula <- function(copula, n) {
  comps <- getAsymExplicitCopulaComps(copula)
  copula1 <- comps$copula1; copula2 <- comps$copula2
  shapes <- comps$shapes
  ## Theorem 2.1, Lemma 2.1, Liebscher (2008, JMA)
  u <- rcopula(copula1, n)
  v <- rcopula(copula2, n)
  d <- copula@dimension
  x <- matrix(NA, n, d)
  for (i in 1:d) {
    x[,i] <- pmax(igfun(u[,i], 1 - shapes[i]), igfun(v[,i], shapes[i]))
  }
  x
}

## setMethod("Afun", signature("asymCopula"), AfunAsymCopula)

setMethod("rcopula", signature("asymExplicitCopula"), rasymExplicitCopula)
setMethod("pcopula", signature("asymExplicitCopula"), pasymExplicitCopula)
setMethod("dcopula", signature("asymExplicitCopula"), dasymExplicitCopula)

