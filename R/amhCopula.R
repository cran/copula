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

setClass("amhCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )


genFunAmh <- function(copula, u) {
  alpha <- copula@parameters[1]
  log((1 - alpha * (1 - u)) / u)
}

genInvAmh <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 - alpha) / (exp(s) - alpha)
}

genFunDer1Amh <- function(copula, u) {
  eval(amhCopula.genfun.expr[1], list(u=u, alpha=copula@parameters[1]))
}
genFunDer2Amh <- function(copula, u) {
  eval(amhCopula.genfun.expr[2], list(u=u, alpha=copula@parameters[1]))
}


amhCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <-   "log((1 - alpha * (1 - u1)) / u1)"
    for (i in 2:n) {
      ui <- paste("u", i, sep="")
      cur <- gsub("u1", ui, expr)
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- gsub("s", expr, "(1 - alpha) / (exp(s) - alpha)")
    parse(text = expr)
  }
  
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  if (dim > 2 && param[1] < 0)
    stop("param can be negative only for dim = 2")
  cdf <- cdfExpr(dim)
  if (dim <= 6)  pdf <- pdfExpr(cdf, dim)
  else pdf <- NULL

  val <- new("amhCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = -1,
             param.upbnd = 1,
             message = "Amh copula family; Archimedean copula")
  val
}


pamhCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  u1 * u2 / (1 - alpha * (1 - u1) * (1 - u2))
}


damhCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  ## bivariate anyway
  u1 <- u[,1]
  u2 <- u[,2]
  (-1 + alpha^2*(-1 + u1 + u2 - u1*u2) - alpha*(-2 + u1 + u2 + u1*u2)) / (-1 + alpha*(-1 + u1)*(-1 + u2))^3 
}


ramhCopula <- function(copula, n) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  val <- matrix(runif(n * dim), nrow = n)
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (val)  ## the limit is independence
  ## Johnson (1987, p.362). Typo V_2 and p?
  ## solve quadratic equation anyway
  u <- runif(n)
  w <- runif(n)
  b <- 1 - u
  A <- w * (alpha * b)^2 - alpha
  B <- (alpha + 1) - 2 * alpha * b * w
  C <- w - 1
  v <- (- B + sqrt(B^2 - 4 * A * C)) / 2 / A
  ## v <- cbind(v, (- B - sqrt(B^2 - 4 * A * C)) / 2 / A) ## this root is not good
  v <- 1 - v
  cbind(u, v)
}

kendallsTauAmhCopula <- function(copula, ...) {
  alpha <- copula@parameters[1]
  ## Nelsen (1999, p.139)
  (3 * alpha - 2) / 3 / alpha - 2 / 3 * (1 - 1/alpha)^2 * log(1 - alpha)
}


setMethod("rcopula", signature("amhCopula"), ramhCopula)
setMethod("pcopula", signature("amhCopula"), pamhCopula)
setMethod("dcopula", signature("amhCopula"), damhCopula)

setMethod("genFun", signature("amhCopula"), genFunAmh)
setMethod("genInv", signature("amhCopula"), genInvAmh)
setMethod("genFunDer1", signature("amhCopula"), genFunDer1Amh)
setMethod("genFunDer2", signature("amhCopula"), genFunDer2Amh)

setMethod("kendallsTau", signature("amhCopula"), kendallsTauAmhCopula)
## setMethod("spearmansRho", signature("amhCopula"), spearmansRhoAmhCopula)
## setMethod("tailIndex", signature("amhCopula"), tailIndexAmhCopula)

## setMethod("calibKendallsTau", signature("amhCopula"), calibKendallsTauAmhCopula)

