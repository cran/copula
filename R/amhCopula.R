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


genFunAmh <- function(copula, u) {
  alpha <- copula@parameters[1]
  log((1 - alpha * (1 - u)) / u)
}

genInvAmh <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 - alpha) / (exp(s) - alpha)
}

genFunDer1Amh <- function(copula, u) {
  eval(amhCopula.genfunDer.expr[1], list(u=u, alpha=copula@parameters[1]))
}
genFunDer2Amh <- function(copula, u) {
  eval(amhCopula.genfunDer.expr[2], list(u=u, alpha=copula@parameters[1]))
}


amhCopula <- function(param, dim = 2L) {
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

##   if (dim > 2 && param[1] < 0)
##     stop("param can be negative only for dim = 2")
  if((dim <- as.integer(dim))> 2) stop("dim can only be 2 for this copula")

  cdf <- cdfExpr(dim)
  pdf <- if (dim <= 6)  pdfExpr(cdf, dim) else NULL

  new("amhCopula",
      dimension = dim,
      parameters = param[1],
      exprdist = c(cdf = cdf, pdf = pdf),
      param.names = "param",
      param.lowbnd = -1,
      param.upbnd = 1,
      message = "Amh copula family; Archimedean copula")
}


pamhCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  u1 * u2 / (1 - alpha * (1 - u1) * (1 - u2))
}


damhCopula <- function(copula, u, log = FALSE, ...) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep.int(if(log) 0 else 1, nrow(u)))
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  ## bivariate anyway
  u1 <- u[,1]
  u2 <- u[,2]
  r <- (-1 + alpha^2*(-1 + u1 + u2 - u1*u2) - alpha*(-2 + u1 + u2 + u1*u2)) /
      (-1 + alpha*(-1 + u1)*(-1 + u2))^3
  if(log) log(r) else r
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

kendallsTauAmhCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## Nelsen (2006, p.172)
  ## range of tau: [(5 - 8 log 2) / 3, 1/3] ~= [-0.1817, 0.3333]
  (3 * alpha - 2) / 3 / alpha - 2 / 3 * (1 - 1/alpha)^2 * log(1 - alpha)
}

## spearmansRhoAmhCopula <- function(copula, ...) {
##   alpha <- copula@parameters[1]
##   ## Nelsen (2006, p.172); need dilog function
##   ## range of rho: 33 - 48 log 2, 4 pi^2 - 39] ~= [-0.2711, 0.4784]
##   12 * (1 + alpha) / alpha^2 * dilog(1 - alpha) - 24 * (1 - alpha) / alpha^2 * log(1 - alpha) - 3 * (alpha + 12) / alpha
## }

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

setMethod("calibKendallsTau", signature("amhCopula"), calibKendallsTauCopula)

