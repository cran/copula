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


AfunGumbel <- function(copula, w) {
  alpha <- copula@parameters[1]
  A <- (w^alpha + (1 - w)^alpha)^(1/alpha)
  ifelse(w == 0 | w == 1, 1, A)
}

AfunDerGumbel <- function(copula, w) {
  alpha <- copula@parameters[1]
  ## deriv(expression((w^alpha + (1 - w)^alpha)^(1/alpha)), "w", hessian=TRUE)
  value <- eval(expression({
    .expr2 <- 1 - w
    .expr4 <- w^alpha + .expr2^alpha
    .expr5 <- 1/alpha
    .expr7 <- .expr5 - 1
    .expr8 <- .expr4^.expr7
    .expr9 <- alpha - 1
    .expr14 <- w^.expr9 * alpha - .expr2^.expr9 * alpha
    .expr15 <- .expr5 * .expr14
    .expr22 <- .expr9 - 1
    .value <- .expr4^.expr5
    .grad <- array(0, c(length(.value), 1L), list(NULL, c("w")))
    .hessian <- array(0, c(length(.value), 1L, 1L), list(NULL,
        c("w"), c("w")))
    .grad[, "w"] <- .expr8 * .expr15
    .hessian[, "w", "w"] <- .expr4^(.expr7 - 1) * (.expr7 * .expr14) *
        .expr15 + .expr8 * (.expr5 * (w^.expr22 * .expr9 * alpha +
        .expr2^.expr22 * .expr9 * alpha))
    attr(.value, "gradient") <- .grad
    attr(.value, "hessian") <- .hessian
    .value
  }), list(alpha=alpha, w=w))
  der1 <- c(attr(value, "gradient"))
  der2 <- c(attr(value, "hessian"))
  data.frame(der1=der1, der2=der2)
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

gumbelCopula <- function(param, dim = 2L) {
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

  cdf <- cdfExpr((dim <- as.integer(dim)))
  pdf <- if (dim <= 6) pdfExpr(cdf, dim) else NULL
  new("gumbelCopula",
             dimension = dim,
             parameters = param[1],
             exprdist = c(cdf = cdf, pdf = pdf),
             param.names = "param",
             param.lowbnd = 1,
             param.upbnd = Inf,
             message = "Gumbel copula family; Archimedean copula; Extreme value copula")
}


rgumbelCopula <- function(copula, n) {
  ## frailty is stable(1,0,0) with 1/alpha
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  ## reduce to indepCopula
  if (alpha - 1 < .Machine$double.eps ^(1/3) ) return(rcopula(indepCopula(dim=dim), n))
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
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  cdf <- copula@exprdist$cdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  val <- eval(cdf)
  pmax(val, 0)
}

dgumbelCopula <- function(copula, u, log=FALSE, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  pdf <- copula@exprdist$pdf
  dim <- copula@dimension
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  if(log) stop("'log=TRUE' not yet implemented")
  eval(pdf)
}

dgumbelCopula.pdf <- function(copula, u) {
  dim <- copula@dimension
  if (dim > 10) stop("Gumbel copula PDF not implemented for dimension > 10.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
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
  if (any(tau < 0)) warning("tau is out of the range [0, 1]")
  ifelse(tau <= 0, 1, 1/(1 - tau))
}


gumbelRhoFun <- function(alpha) {
  ss <- .gumbelRho$ss
  forwardTransf <- .gumbelRho$trFuns$forwardTransf
  valFun <- .gumbelRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(valFun(theta))
}

gumbelRhoDer <- function(alpha) {
  ss <- .gumbelRho$ss
  forwardTransf <- .gumbelRho$trFuns$forwardTransf
  forwardDer <- .gumbelRho$trFuns$forwardDer
  valFun <- .gumbelRho$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  c(valFun(theta, 1)) * forwardDer(alpha, ss)
}

spearmansRhoGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  gumbelRhoFun(alpha)
}


calibSpearmansRhoGumbelCopula <- function(copula, rho) {
  if (any(rho < 0)) warning("rho is out of the range [0, 1]")
  gumbelRhoInv <- approxfun(x = .gumbelRho$assoMeasFun$fm$ysmth,
                            y = .gumbelRho$assoMeasFun$fm$x, rule = 2)

  ss <- .gumbelRho$ss
  theta <- gumbelRhoInv(rho)
  .gumbelRho$trFuns$backwardTransf(theta, ss)
}


tauDerGumbelCopula <- function(copula) {
  return( 1 / copula@parameters^2 )
}

rhoDerGumbelCopula <- function(copula) {
  alpha <- copula@parameters[1]
  gumbelRhoDer(alpha)
}

setMethod("rcopula", signature("gumbelCopula"), rgumbelCopula)
setMethod("pcopula", signature("gumbelCopula"), pgumbelCopula)
setMethod("dcopula", signature("gumbelCopula"), dgumbelCopula.pdf)

setMethod("Afun", signature("gumbelCopula"), AfunGumbel)
setMethod("AfunDer", signature("gumbelCopula"), AfunDerGumbel)

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
