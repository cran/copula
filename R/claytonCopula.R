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


iPsiClayton <- function(copula, u) {
  alpha <- copula@parameters[1]
  ## (u^(-alpha) - 1) / alpha ## This is in Nelson Table 4.1; alpha in denom redundant
  u^(-alpha) - 1
}

psiClayton <- function(copula, s) {
  alpha <- copula@parameters[1]
  ## (1 + alpha * s)^(-1/alpha) ## corresponding to the comment above
  (1 + s)^(-1/alpha)
}

claytonCopula <- function(param = NA_real_, dim = 2L) {
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

  if ((dim <- as.integer(dim)) > 2 && param[1] < 0)
    stop("param can be negative only for dim = 2")
  cdf <- cdfExpr(dim)
  pdf <- if (dim <= 6) pdfExpr(cdf, dim) # else NULL
  new("claytonCopula",
      dimension = dim,
      parameters = param[1],
      exprdist = c(cdf = cdf, pdf = pdf),
      param.names = "param",
      param.lowbnd = if(dim == 2) -1 else 0,
      param.upbnd = Inf,
      fullname = "Clayton copula family; Archimedean copula")
}

rclaytonBivCopula <- function(n, copula) {
  val <- cbind(runif(n), runif(n))
  alpha <- copula@parameters[1]
  ## This implementation is confirmed by Splus module finmetrics
  val[,2] <- (val[,1]^(-alpha) * (val[,2]^(-alpha/(alpha + 1)) - 1) + 1)^(-1/alpha)
  ## Frees and Valdez (1998, p.11): wrong! Can be checked by sample Kendall' tau.
  ## Their general expression for $F_k$ is all right for $k > 2$.
  ## But for dimension 2, they are missing $\phi'$ in the numerator and
  ## the denominator should be 1. So their formula for $U_2$ on p.11 is incorrect.
  val
}


rclaytonCopula <- function(n, copula) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha - 0) < .Machine$double.eps ^ (1/3))
    return(rCopula(n, indepCopula(dim)))
  if (dim == 2) return (rclaytonBivCopula(n, copula))
  ## gamma frailty
  val <- matrix(runif(n * dim), nrow = n)
  if (abs(alpha) <= 100 * .Machine$double.eps)
    return (val)  ## the limit is independence
  ## gam <- rgamma(n, shape = 1/alpha, rate = 1/alpha) ## fixed from rate = 1
  gam <- rgamma(n, shape = 1/alpha, rate = 1) ## fixed from rate = 1
  gam <- matrix(gam, nrow = n, ncol = dim)
  psi(copula, - log(val) / gam)
}


pclaytonCopula <- function(copula, u) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  cdf <- copula@exprdist$cdf
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  pmax(eval(cdf), 0)
}


## Nowhere used (!)
dclaytonCopula <- function(copula, u, ...) {
  dim <- copula@dimension
  if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (rep(1, nrow(u)))
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  val <- c(eval(copula@exprdist$pdf))
  val[apply(u, 1, function(v) any(v < 0))] <- 0
  val[apply(u, 1, function(v) any(v > 1))] <- 0
##   if (alpha < 0) {
##     cdf <- pCopula(u, copula)
##     bad <- cdf == 0
##     val[bad] <- 0
##   }
  if(log) log(val) else val
}

dclaytonCopula.pdf <- function(u, copula, log=FALSE) {
  dim <- copula@dimension
  if (dim > 10) stop("Clayton copula PDF not implemented for dimension > 10.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9)
    return(rep.int(if(log) 0 else 1, nrow(u)))
  val <- pmax(c(eval(claytonCopula.pdf.algr[dim])),0)
  ## clean up
  val[apply(u, 1, function(v) any(v < 0))] <- 0
  val[apply(u, 1, function(v) any(v > 1))] <- 0

  val[apply(u, 1, function(v) any(v == 0) & any(v > 0))] <- 0

  ## if (alpha > 0)
  ## else
  if(log) log(val) else val
}


tauClaytonCopula <- function(copula) {
  alpha <- copula@parameters[1]
  alpha / (alpha + 2)
}

iTauClaytonCopula <- function(copula, tau) {
  2 * tau / (1 - tau)
}

claytonRhoFun <- function(alpha) {
  ss <- .claytonRhoNeg$ss
  forwardTransf <- .claytonRhoNeg$trFuns$forwardTransf
  valFunNeg <- .claytonRhoNeg$assoMeasFun$valFun
  valFunPos <- .claytonRhoPos$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  as.vector(if(alpha <= 0) valFunNeg(theta) else valFunPos(theta))
}

claytonRhoDer <- function(alpha) {
  ss <- .claytonRhoNeg$ss
  forwardTransf <- .claytonRhoNeg$trFuns$forwardTransf
  forwardDer <- .claytonRhoNeg$trFuns$forwardDer
  valFunNeg <- .claytonRhoNeg$assoMeasFun$valFun
  valFunPos <- .claytonRhoPos$assoMeasFun$valFun
  theta <- forwardTransf(alpha, ss)

  as.vector(if(alpha <= 0) valFunNeg(theta, 1) else valFunPos(theta, 1)) * forwardDer(alpha, ss)
}

rhoClaytonCopula <- function(copula) {
  alpha <- copula@parameters[1]
  claytonRhoFun(alpha)
}


iRhoClaytonCopula <- function(copula, rho) {
  claytonRhoInvNeg <- approxfun(x = .claytonRhoNeg$assoMeasFun$fm$ysmth,
                                y = .claytonRhoNeg$assoMeasFun$fm$x)

  claytonRhoInvPos <- approxfun(x = .claytonRhoPos$assoMeasFun$fm$ysmth,
                                y = .claytonRhoPos$assoMeasFun$fm$x)

  ss <- .claytonRhoNeg$ss
  theta <- if(rho <= 0) claytonRhoInvNeg(rho) else claytonRhoInvPos(rho)
  .claytonRhoPos$trFuns$backwardTransf(theta, ss)
}

dRhoClaytonCopula <- function(copula) {
  alpha <- copula@parameters[1]
  claytonRhoDer(alpha)
}

tailIndexClaytonCopula <- function(copula) {
  alpha <- copula@parameters
  c(lower= if(alpha > 0) 2 ^ (-1/alpha) else 0,
    upper= 0)
}

dTauClaytonCopula <- function(copula) {
  return( 2 / (copula@parameters+2)^2 )
}

pMatClayton <- function (u, copula, ...) {
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && !copClayton@paraConstr(th)) # for now, .. to support negative tau
        pclaytonCopula(copula, u=u)
    else
        pacopula(u, copClayton, theta=copula@parameters, ...)
}

dMatClayton <- function (u, copula, log = FALSE, ...) {
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && !copClayton@paraConstr(th)) # for now, .. to support negative tau
        dclaytonCopula.pdf(u, copula, log=log)
    else
        copClayton@dacopula(u, theta=copula@parameters, log=log, ...)
}

setMethod("rCopula", signature("numeric", "claytonCopula"), rclaytonCopula)

setMethod("pCopula", signature("numeric", "claytonCopula"),
	  function (u, copula, ...)
	  ## was  pclaytonCopula
          pMatClayton(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("pCopula", signature("matrix", "claytonCopula"), pMatClayton)

setMethod("dCopula", signature("numeric", "claytonCopula"),
	  function (u, copula, ...)
	  ## was  dclaytonCopula.pdf
          dMatClayton(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("dCopula", signature("matrix", "claytonCopula"), dMatClayton)


setMethod("iPsi", signature("claytonCopula"), iPsiClayton)
## FIXME {negative tau}
## setMethod("iPsi", signature("claytonCopula"),
## 	  function(copula, u) copClayton@iPsi(u, theta=copula@parameters))
setMethod("psi", signature("claytonCopula"), psiClayton)
## FIXME {negative tau}
## setMethod("psi", signature("claytonCopula"),
## 	  function(copula, s) copClayton@psi(t=s, theta=copula@parameters))

setMethod("diPsi", signature("claytonCopula"),
	  function(copula, u, degree=1, log=FALSE, ...)
      {
	  s <- if(log || degree %% 2 == 0) 1. else -1.
	  s* copClayton@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
      })



setMethod("tau", signature("claytonCopula"), tauClaytonCopula)
setMethod("rho", signature("claytonCopula"), rhoClaytonCopula)
setMethod("tailIndex", signature("claytonCopula"), tailIndexClaytonCopula)

setMethod("iTau", signature("claytonCopula"), iTauClaytonCopula)
setMethod("iRho", signature("claytonCopula"), iRhoClaytonCopula)

setMethod("dTau", signature("claytonCopula"), dTauClaytonCopula)
setMethod("dRho", signature("claytonCopula"), dRhoClaytonCopula)
