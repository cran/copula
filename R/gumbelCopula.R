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


AfunGumbel <- function(copula, w) {
  alpha <- copula@parameters[1]
  (w^alpha + (1 - w)^alpha)^(1/alpha)
}

genFunGumbel <- function(copula, u) {
  alpha <- copula@parameters[1]
  ( - log(u))^alpha
}

genInvGumbel <- function(copula, s) {
  alpha <- copula@parameters[1]
  exp( -s^(1 / alpha) )
}

genFunDer1Gumbel <- function(copula, u) {
  eval(gumbelCopula.genfunDer.expr[1], list(u=u, alpha=copula@parameters[1]))
}


genFunDer2Gumbel <- function(copula, u) {
  eval(gumbelCopula.genfunDer.expr[2], list(u=u, alpha=copula@parameters[1]))
}

gumbelCopula <- function(param, dim = 2) {
  ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    expr <- "( - log(u1))^alpha"
    for (i in 2:n) {
      cur <- paste( "(-log(u", i, "))^alpha", sep = "")
      expr <- paste(expr, cur, sep=" + ")
    }
    expr <- paste("exp(- (", expr, ")^ (1/alpha))")
    parse(text = expr)
  }
  
  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }

  cdf <- cdfExpr(dim)
  if (dim <= 6)  pdf <- pdfExpr(cdf, dim)
  else pdf <- NULL
  val <- new("gumbelCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = 1,
             param.upbnd = Inf,
             message = "Gumbel copula family; Archimedean copula; Extreme value copula")
  val
}

## rgumbelBivCopula <- function(copula, n) {
##   ## Taken from finmetrics for comparison purpose, but I don't understand it
##   ## because it is not well-documented.
##   ## generate z from distribution h
##   ## using rejection method
##   zrand <- vector(length = n, mode = "numeric")
##   xgenerated <- vector(length = n, mode = "logical")
##   lg <- n
##   ##set c to 1/2 true for symmetric copula
##   cc <- Hderiv(copula, 1/2)
##   stopnow <- F
##   while(!stopnow) {
##     usamp <- runif(lg)
##     ysamp <- runif(lg)
##     ykeep <- (usamp <= Hderiv(copula, ysamp)/cc)
##     ##indexes of xrand vector to store "good" generated values  
##     toind <- (c(1:n)[!xgenerated])[ykeep]
##     zrand[toind] <- ysamp[ykeep]
##     xgenerated[toind] <- T
##     ##number of rvs left to generate:
##     lg <- n - sum(xgenerated)
##     if(lg == 0)
##       stopnow <- T
##   }
##   z <- zrand
##   u <- runif(n)
##   ## cat(length(z), length(AsecondDer(copula,z)))
##   pz <- (z * (1 - z) * AsecondDer(copula, z))/Hderiv(copula, z)/A(copula, z)
##   w <- NULL
##   w[1:n] <- 0
##   nn <- sum(u <= pz)
##   w[u <= pz] <- runif(nn)
##   w[u > pz] <- runif(n - nn) * runif(n - nn)
##   xx <- exp((z * log(w))/A(copula, z))
##   yy <- exp(((1 - z) * log(w))/A(copula, z))
##   val <- list(x = xx, y = yy)
##   val
## } 

rgumbelCopula <- function(copula, n) {
  ## frailty is stable(1,0,0) with 1/alpha
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  b <- 1/alpha
  ## stable (b, 1), 0 < b < 1, Chambers, Mallows, and Stuck 1976, JASA, p.341
  fr <- rPosStable(n, b)
  fr <- matrix(fr, nrow=n, ncol=dim)
  ## now gumbel copula
  val <- matrix(runif(dim * n), nrow = n)
  genInv(copula, - log(val) / fr)
}


pgumbelCopula <- function(copula, u) {
  dim <- copula@dimension
  if (is.vector(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  val <- eval(cdf)
  pmax(val, 0)
}

dgumbelCopula <- function(copula, u) {
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  eval(pdf)
}

dgumbelCopula.pdf <- function(copula, u) {
  dim <- copula@dimension
  if (dim > 10) stop("Gumbel copula PDF not implemented for dimension > 10.")
  if (is.vector(u)) u <- matrix(u, nrow = 1)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  c(eval(gumbelCopula.pdf.algr[dim]))
}

kendallsTauGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  1 - 1/alpha
}

tailIndexGumbelCopula <- function(copula, ...) {
  alpha <- copula@parameters
  upper <- 2 - 2^(1/alpha)
  c(lower=0, upper=upper)
}


calibKendallsTauGumbelCopula <- function(copula, tau) {
  1/(1 - tau)
}

spearmansRhoGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha > 15) spearmansRhoCopula(copula)
  else spearmansRhoGumbelCopula.tr(alpha)
}

calibSpearmansRhoGumbelCopula <- function(copula, rho) {
  if (rho > 0.96366) calibSpearmansRhoCopula(copula, rho)
  else calibSpearmansRhoGumbelCopula.tr(rho)
}

tauDerGumbelCopula <- function(copula)
  {
    return( 1 / copula@parameters^2 )
  }

rhoDerGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha > 15) stop("not implemented yet.")
  else spearmansRhoDerGumbelCopula.tr(alpha)
}

setMethod("rcopula", signature("gumbelCopula"), rgumbelCopula)
setMethod("pcopula", signature("gumbelCopula"), pgumbelCopula)
setMethod("dcopula", signature("gumbelCopula"), dgumbelCopula.pdf)

setMethod("Afun", signature("gumbelCopula"), AfunGumbel)

setMethod("genFun", signature("gumbelCopula"), genFunGumbel)
setMethod("genInv", signature("gumbelCopula"), genInvGumbel)

setMethod("genFunDer1", signature("gumbelCopula"), genFunDer1Gumbel)
setMethod("genFunDer2", signature("gumbelCopula"), genFunDer2Gumbel)

setMethod("kendallsTau", signature("gumbelCopula"), kendallsTauGumbelCopula)
setMethod("spearmansRho", signature("gumbelCopula"), spearmansRhoGumbelCopula)
setMethod("tailIndex", signature("gumbelCopula"), tailIndexGumbelCopula)

setMethod("calibKendallsTau", signature("gumbelCopula"), calibKendallsTauGumbelCopula)
setMethod("calibSpearmansRho", signature("gumbelCopula"), calibSpearmansRhoGumbelCopula)

setMethod("rhoDer", signature("gumbelCopula"), rhoDerGumbelCopula)
setMethod("tauDer", signature("gumbelCopula"), tauDerGumbelCopula)
