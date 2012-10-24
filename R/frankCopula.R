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


iPsiFrank <- function(copula, u) {
  alpha <- copula@parameters[1]
  - log( (exp(- alpha * u) - 1) / (exp(- alpha) - 1))
}

psiFrank <- function(copula, s) {
  alpha <- copula@parameters[1]
  -1/alpha * log(1 + exp(-s) * (exp(-alpha) - 1))
}

## psiDerFrank <- function(copula, s, n) {
##   eval(psiDerFrank.expr[n + 1], list(s=s, alpha=copula@parameters[1]))
## }

frankCopula <- function(param = NA_real_, dim = 2L) {
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
  new("frankCopula",
      dimension = dim,
      parameters = param[1],
      exprdist = c(cdf = cdf, pdf = pdf),
      param.names = "param",
      param.lowbnd = -Inf,
      param.upbnd = Inf,
      fullname = "Frank copula family; Archimedean copula")
}

rfrankBivCopula <- function(n, copula) {
  val <- cbind(runif(n), runif(n))
  ## to fix numerical rounding problems for alpha >35 but not for alpha < -35
  alpha <- - abs(copula@parameters[1])
  val[,2] <- -1/alpha * log(1 + val[,2] * (1 - exp(-alpha)) / (exp(-alpha * val[,1]) * (val[,2] - 1) - val[,2])) ## reference: Joe (1997, p.147)
  if (copula@parameters[1] > 0) val[,2] <- 1 - val[,2]
  val
}

## rfrankBivCopula <- function(n, copula) {
##   delta <- copula@parameters[1]
##   q <- runif(n); u <- runif(n)
##   v <- - 1/ delta * log( 1 - (1 - exp(-delta)) / ((1 / q - 1) * exp(- delta * u) + 1))
##   cbind(u, v)
## }

rfrankCopula <- function(n, copula) {
  dim <- copula@dimension
  alpha <- copula@parameters[1]
  if (abs(alpha - 0) < .Machine$double.eps ^ (1/3))
    return(rCopula(n, indepCopula(dim)))
##   if (abs(alpha) <= .Machine$double.eps^.9)
##     return (matrix(runif(n * dim), nrow = n))
  if (dim == 2) return (rfrankBivCopula(n, copula))
  ## the frailty is a log series distribution with a = 1 - exp(-alpha)
  fr <- rlogseries(n, 1 - exp(-alpha))
  fr <- matrix(fr, nrow = n, ncol = dim)
  val <- matrix(runif(dim * n), nrow = n)
  psi(copula, - log(val) / fr)
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

dfrankCopula <- function(u, copula, log=FALSE, ...) {
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
##   s <- apply(iPsiFrank(copula, u), 1, sum)
##   pdf <- psiDerFrank(copula, s, copula@dimension) *
##     apply(iPsiDerFrank(copula, u, 1), 1, prod)
##   pdf
## }

dfrankCopula.pdf <- function(u, copula, log=FALSE) {
  dim <- copula@dimension
  if (dim > 6) stop("Frank copula PDF not implemented for dimension > 6.")
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
  alpha <- copula@parameters[1]
  ## FIXME: improve log-case
  if(log)
    log(c(eval(frankCopula.pdf.algr[dim])))
  else  c(eval(frankCopula.pdf.algr[dim]))
}


tauFrankCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha == 0) return (0)
  1 - 4 / alpha * (1 - debye1(alpha))
}

rhoFrankCopula <- function(copula) {
  alpha <- copula@parameters[1]
  if (alpha == 0) 0
  else
      ## Genest (1987), "Frank's family of bivariate distributions" (Biometrika, 7, 549--555)
      1 - 12/alpha * (debye1(alpha) - debye2(alpha))
}

dTauFrankCopula <- function(copula) {
  alpha <- copula@parameters
  4/alpha^2 + 4/(alpha * expm1(alpha)) - 8/alpha^2 * debye1(alpha)
}

dRhoFrankCopula <- function(copula) {
  alpha <- copula@parameters
  12 / (alpha * expm1(alpha)) + (-36 * debye2(alpha) + 24 * debye1(alpha))/ alpha^2
}


pMatFrank <- function (u, copula, ...) {
    ## was  pfrankCopula
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && !copFrank@paraConstr(th)) # for now, .. to support negative tau
        pfrankCopula(copula, u=u)
    else
        pacopula(u, copFrank, theta=copula@parameters, ...)
}

dMatFrank <- function (u, copula, log = FALSE, ...) {
    ## was  dfrankCopula.pdf
    stopifnot(ncol(u) == (d <- copula@dimension))
    th <- copula@parameters
    if(d == 2 && th < 0) # for now, copFrank does not yet support negative tau
        dfrankCopula.pdf(u, copula, log=log)
    else
        copFrank@dacopula(u, theta=th, log=log, ...)
}
setMethod("rCopula", signature("numeric", "frankCopula"), rfrankCopula)

setMethod("pCopula", signature("numeric", "frankCopula"),
	  function (u, copula, ...)
	  pMatFrank(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("pCopula", signature("matrix", "frankCopula"), pMatFrank)

setMethod("dCopula", signature("numeric", "frankCopula"),
	  function (u, copula, ...)
	  dMatFrank(matrix(u, ncol = dim(copula)), copula, ...))
setMethod("dCopula", signature("matrix", "frankCopula"), dMatFrank)

setMethod("iPsi", signature("frankCopula"), iPsiFrank)
## FIXME {negative tau}
## setMethod("iPsi", signature("frankCopula"),
## 	  function(copula, u) copFrank@iPsi(u, theta=copula@parameters))
setMethod("psi", signature("frankCopula"), psiFrank)
## FIXME {negative tau}
## setMethod("psi", signature("frankCopula"),
## 	  function(copula, s) copFrank@psi(t=s, theta=copula@parameters))

## setMethod("psiDer", signature("frankCopula"), psiDerFrank)
setMethod("diPsi", signature("frankCopula"),
	  function(copula, u, degree=1, log=FALSE, ...)
      {
	  s <- if(log || degree %% 2 == 0) 1. else -1.
	  s* copFrank@absdiPsi(u, theta=copula@parameters, degree=degree, log=log, ...)
      })



setMethod("tau", signature("frankCopula"), tauFrankCopula)
setMethod("rho", signature("frankCopula"), rhoFrankCopula)
setMethod("tailIndex", signature("frankCopula"), function(copula) c(lower=0, upper=0))

setMethod("iTau", signature("frankCopula"),
	  function(copula, tau, tol = 1e-7) copFrank@iTau(tau, tol=tol))
setMethod("iRho", signature("frankCopula"), iRhoCopula)

setMethod("dRho", signature("frankCopula"), dRhoFrankCopula)
setMethod("dTau", signature("frankCopula"), dTauFrankCopula)
