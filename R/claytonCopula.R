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

setClass("claytonCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )


genFunClayton <- function(copula, u) {
  alpha <- copula@parameters[1]
  (u^(-alpha) - 1) / alpha
}

genInvClayton <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 + alpha * s)^(-1/alpha) 
}

genFunDer1Clayton <- function(copula, u) {
  alpha <- copula@parameters[1]
  eval(claytonCopula.genfun.expr[1])
}

genFunDer2Clayton <- function(copula, u) {
  alpha <- copula@parameters[1]
  eval(claytonCopula.genfun.expr[2])
}


claytonCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <- "u1^(-alpha) - 1"
    for (i in 2:n) {
      cur <- paste( "u", i, "^(-alpha) - 1", sep = "")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("(1 + (", expr, "))^ (-1/alpha)")
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
  val <- new("claytonCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = ifelse(dim == 2, -1, 0),
             param.upbnd = Inf,
             message = "Clayton copula family; Archimedean copula")
  val
}

rclaytonBivCopula <- function(copula, n) {
  val <- cbind(runif(n), runif(n))
  alpha <- copula@parameters[1]
  val[,2] <- (val[,1]^(-alpha) * (val[,2]^(-alpha/(alpha + 1)) - 1) + 1)^(-1/alpha)
  val
}


rclaytonCopula <- function(copula, n) {
  dim <- copula@dimension
  if (dim == 2) return (rclaytonBivCopula(copula, n))
  ## gamma frailty
  alpha <- copula@parameters[1]
  val <- matrix(runif(n * dim), nrow = n)
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (val)  ## the limit is independence
  gam <- rgamma(n, shape = 1/alpha , rate = 1)
  gam <- matrix(gam, nrow = n, ncol = dim)
  genInv(copula, - log(val) / gam)
}


pclaytonCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  cdf <- copula@exprdist$cdf
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  val <- eval(cdf)
  pmax(val, 0)
}


dclaytonCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  pdf <- copula@exprdist$pdf
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  val <- eval(pdf)
  val[apply(u, 1, function(v) any(v < 0))] <- 0
  val[apply(u, 1, function(v) any(v > 1))] <- 0
##   if (alpha < 0) {
##     cdf <- pcopula(copula, u)
##     bad <- cdf == 0
##     val[bad] <- 0
##   }
  val
}

dclaytonCopula.pdf <- function(copula, u) {
  dim <- copula@dimension
  if (dim > 10) stop("Clayton copula PDF not implemented for dimension > 10.")
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  val <- eval(claytonCopula.pdf.algr[dim])
  ## clean up
  val[apply(u, 1, function(v) any(v < 0))] <- 0
  val[apply(u, 1, function(v) any(v > 1))] <- 0
  
  val[apply(u, 1, function(v) any(v == 0) & any(v > 0))] <- 0
  
  ## if (alpha > 0)
  ## else 
  val
}


kendallsTauClaytonCopula <- function(copula) {
  alpha <- copula@parameters[1]
  alpha / (alpha + 2)
}

calibKendallsTauClaytonCopula <- function(copula, tau) {
  2 * tau / (1 - tau)
}

tailIndexClaytonCopula <- function(copula, ...) {
  upper <- 0
  alpha <- copula@parameters
  lower <- ifelse(alpha > 0, 2 ^ (-1/alpha), 0)
  c(lower=lower, upper=upper)
}


setMethod("rcopula", signature("claytonCopula"), rclaytonCopula)
setMethod("pcopula", signature("claytonCopula"), pclaytonCopula)
setMethod("dcopula", signature("claytonCopula"), dclaytonCopula)

setMethod("genFun", signature("claytonCopula"), genFunClayton)
setMethod("genInv", signature("claytonCopula"), genInvClayton)
setMethod("genFunDer1", signature("claytonCopula"), genFunDer1Clayton)
setMethod("genFunDer2", signature("claytonCopula"), genFunDer2Clayton)

setMethod("kendallsTau", signature("claytonCopula"), kendallsTauClaytonCopula)
#setMethod("spearmansRho", signature("claytonCopula"), spearmansRhoCopula)
setMethod("tailIndex", signature("claytonCopula"), tailIndexClaytonCopula)

setMethod("calibKendallsTau", signature("claytonCopula"), calibKendallsTauClaytonCopula)

