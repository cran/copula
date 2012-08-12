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


iPsiAmh <- function(copula, u) {
  alpha <- copula@parameters[1]
  log((1 - alpha * (1 - u)) / u)
}

psiAmh <- function(copula, s) {
  alpha <- copula@parameters[1]
  (1 - alpha) / (exp(s) - alpha)
}

amhCopula <- function(param = NA_real_, dim = 2L) {
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
      param.lowbnd = -1,# 0 for tau >= 0
      param.upbnd = 1,
      fullname = "Amh copula family; Archimedean copula")
}


pamhCopula <- function(copula, u) {
  dim <- copula@dimension
  ## if(!is.matrix(u)) u <- matrix(u, ncol = dim)
  alpha <- copula@parameters[1]
  if (abs(alpha) <= .Machine$double.eps^.9) return (apply(u, 1, prod))
  ## for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  u1 <- u[,1]
  u2 <- u[,2]
  u1 * u2 / (1 - alpha * (1 - u1) * (1 - u2))
}


damhCopula <- function(u, copula, log = FALSE, ...) {
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


ramhCopula <- function(n, copula) {
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

iTauAmhCopula <- function(copula, tau) {
  tauMin <- (5 - 8*log(2)) / 3
  iSml <- tau < tauMin
  iLrg <- tau > 1/3
  if(any(out <- iSml | iLrg))
      warning("tau is out of the range [(5 - 8 log 2) / 3, 1/3] ~= [-0.1817, 0.3333]")
  ## A _faster_ version of  r <- ifelse(tau < tauMin, -1, ifelse(tau > 1/3, 1, ....)):
  r <- tau
  r[iSml] <- -1.
  r[iLrg] <- +1.
  if(any(in. <- !out))
      r[in.] <- iTauCopula(copula, tau[in.])
  r
}

rhoAmhCopula <- function(copula) {
  alpha <- copula@parameters[1]
  ## Nelsen (2006, p.172); need dilog function, where his dilog(x) = Li_2(1-x) = polylog(1-x, 2)
  ## range of rho: 33 - 48 log 2, 4 pi^2 - 39] ~= [-0.2711, 0.4784]
  12 * (1 + alpha) / alpha^2 * polylog(alpha,2) - 24 * (1 - alpha) / alpha^2 * log1p(- alpha) -
      3 * (alpha + 12) / alpha
}

pMatAmh <- function (u, copula, ...) {
    ## was pamhCopula
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && !copAMH@paraConstr(th)) # for now, .. to support negative tau
        pamhCopula(copula, u)
    else
        pacopula(u, copAMH, theta=copula@parameters, ...)
}

dMatAmh <- function (u, copula, log = FALSE, ...) {
    ## was  damhCopula
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && !copAMH@paraConstr(th)) # for now, .. to support negative tau
        damhCopula(u, copula, log=log)
    else
        copAMH@dacopula(u, theta=copula@parameters, log=log, ...)
}

setMethod("rCopula", signature("numeric", "amhCopula"), ramhCopula)

setMethod("pCopula", signature("numeric", "amhCopula"),
	  function (u, copula, ...)
          pMatAmh(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("pCopula", signature("matrix", "amhCopula"), pMatAmh)

setMethod("dCopula", signature("numeric", "amhCopula"),
	  function (u, copula, ...)
          dMatAmh(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("dCopula", signature("matrix", "amhCopula"), dMatAmh)

setMethod("iPsi", signature("amhCopula"), iPsiAmh)
## setMethod("iPsi", signature("amhCopula"),
## 	  function(copula, u) copAMH@iPsi(u, theta=copula@parameters))
setMethod("psi", signature("amhCopula"), psiAmh)
## FIXME
## setMethod("psi", signature("amhCopula"),
## 	  function(copula, s) copAMH@psi(t=s, theta=copula@parameters))
setMethod("diPsi", signature("amhCopula"),
	  function(copula, u, degree=1, log=FALSE, ...) {
              s <- if(log || degree %% 2 == 0) 1. else -1.
              s* copAMH@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
      })


setMethod("tau", signature("amhCopula"), function(copula) tauAMH(copula@parameters[1]))

setMethod("rho", signature("amhCopula"), rhoAmhCopula)
setMethod("tailIndex", signature("amhCopula"), function(copula) c(lower=0, upper=0))

setMethod("iTau", signature("amhCopula"), iTauAmhCopula)

