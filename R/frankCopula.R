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


genFunFrank <- function(copula, u) {
  alpha <- copula@parameters[1]
  - log( (exp(- alpha * u) - 1) / (exp(- alpha) - 1))
}

genInvFrank <- function(copula, s) {
  alpha <- copula@parameters[1]
  -1/alpha * log(1 + exp(-s) * (exp(-alpha) - 1))
}

genFunDer1Frank <- function(copula, u) {
  eval(frankCopula.genfunDer.expr[1], list(u=u, alpha=copula@parameters[1]))
}

genFunDer2Frank <- function(copula, u) {
  eval(frankCopula.genfunDer.expr[2], list(u=u, alpha=copula@parameters[1]))
}

## genInvDerFrank <- function(copula, s, n) {
##   eval(genInvDerFrank.expr[n + 1], list(s=s, alpha=copula@parameters[1]))
## }

frankCopula <- function(param, dim = 2L) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <-   "- log( (exp(- alpha * u1) - 1) / (exp(- alpha) - 1) )"
    for (i in 2:n) {
      cur <- paste("- log( (exp(- alpha * u", i, ") - 1) / (exp(- alpha) - 1))", sep="")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("-1/alpha * log(1 + exp(-(", expr, ")) * (exp(-alpha) - 1))")
    parse(text = expr)
  }

  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  if ((dim <- as.integer(dim)) > 2 && param[1] < 0)
    stop("param can be negative only for dim = 2")
  cdf <- cdfExpr(dim)
  if (dim <= 6)  pdf <- pdfExpr(cdf, dim)
  else pdf <- NULL
  val <- new("frankCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = -Inf,
             param.upbnd = Inf,
             message = "Frank copula family; Archimedean copula")
  val
}

rfrankBivCopula <- function(copula, n) {
  val <- cbind(runif(n), runif(n))
  ## to fix numerical rounding problems for alpha >35 but not for alpha < -35
  alpha <- - abs(copula@parameters[1])
  val[,2] <- -1/alpha * log(1 + val[,2] * (1 - exp(-alpha)) / (exp(-alpha * val[,1]) * (val[,2] - 1) - val[,2])) ## reference: Joe (1997, p.147)
  if (copula@parameters[1] > 0) val[,2] <- 1 - val[,2]
  val
}

## rfrankBivCopula <- function(copula, n) {
##   delta <- copula@parameters[1]
##   q <- runif(n); u <- runif(n)
##   v <- - 1/ delta * log( 1 - (1 - exp(-delta)) / ((1 / q - 1) * exp(- delta * u) + 1))
##   cbind(u, v)
## }

rfrankCopula <- function(copula, n) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha - 0) < .Machine$double.eps ^ (1/3))
    return(rcopula(indepCopula(dim), n))
##   if (abs(alpha) <= .Machine$double.eps^.9)
##     return (matrix(runif(n * dim), nrow = n))
  if (dim == 2) return (rfrankBivCopula(copula, n))
  ## the frailty is a log series distribution with a = 1 - exp(-alpha)
  fr <- rlogseries(n, 1 - exp(-alpha))
  fr <- matrix(fr, nrow = n, ncol = dim)
  val <- matrix(runif(dim * n), nrow = n)
  genInv(copula, - log(val) / fr)
}


pfrankCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  eval(cdf)
}

dfrankCopula <- function(copula, u, log=FALSE, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  if(log) stop("'log=TRUE' not yet implemented")
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  val <- eval(pdf)
#  val[apply(u, 1, function(v) any(v <= 0))] <- 0
#  val[apply(u, 1, function(v) any(v >= 1))] <- 0
  val
}

## dfrankCopula.expr <- function(copula, u) {
##   if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
##   s <- apply(genFunFrank(copula, u), 1, sum)
##   pdf <- genInvDerFrank(copula, s, copula@dimension) *
##     apply(genFunDerFrank(copula, u, 1), 1, prod)
##   pdf
## }

dfrankCopula.pdf <- function(copula, u) {
  dim <- copula@dimension
  if (dim > 6) stop("Frank copula PDF not implemented for dimension > 6.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  c(eval(frankCopula.pdf.algr[dim]))
}

kendallsTauFrankCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  1 - 4 / alpha * (1 - debye1(alpha))
}

spearmansRhoFrankCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
    1 - 12/alpha * (debye1(alpha) - debye2(alpha))
}

tailIndexFrankCopula <- function(copula) {
  c(lower=0, upper=0)
}

tauDerFrankCopula <- function(copula) {
  alpha <- copula@parameters
  return( 4/alpha^2 + 4/(alpha * (exp(alpha) - 1)) - 8/alpha^2 * debye1(alpha) )
}

rhoDerFrankCopula <- function(copula) {
  alpha <- copula@parameters
  return( 12 / (alpha * (exp(alpha) - 1)) - 36 / alpha^2 * debye2(alpha) + 24 / alpha^2 * debye1(alpha) )
}


setMethod("rcopula", signature("frankCopula"), rfrankCopula)
setMethod("pcopula", signature("frankCopula"), pfrankCopula)
setMethod("dcopula", signature("frankCopula"), dfrankCopula.pdf)

setMethod("genFun", signature("frankCopula"), genFunFrank)
setMethod("genInv", signature("frankCopula"), genInvFrank)
setMethod("genFunDer1", signature("frankCopula"), genFunDer1Frank)
setMethod("genFunDer2", signature("frankCopula"), genFunDer2Frank)
## setMethod("genInvDer", signature("frankCopula"), genInvDerFrank)


setMethod("kendallsTau", signature("frankCopula"), kendallsTauFrankCopula)
setMethod("spearmansRho", signature("frankCopula"), spearmansRhoFrankCopula)
setMethod("tailIndex", signature("frankCopula"), tailIndexFrankCopula)

setMethod("calibKendallsTau", signature("frankCopula"), calibKendallsTauCopula)
setMethod("calibSpearmansRho", signature("frankCopula"), calibSpearmansRhoCopula)

setMethod("rhoDer", signature("frankCopula"), rhoDerFrankCopula)
setMethod("tauDer", signature("frankCopula"), tauDerFrankCopula)
