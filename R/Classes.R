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



# validCopula <- function(object) {
#   dim <- object@dimension
#   if (dim != as.integer(dim))
#     return("dim must be integer")
#   if (dim < 2)
#     return("dim must be >= 2")
#   param <- object@parameters
#   upper <- object@param.upbnd
#   lower <- object@param.lowbnd
#   if (length(param) != length(upper))
#     return("Parameter and upper bound have non-equal length")
#   if (length(param) != length(lower))
#     return("Parameter and lower bound have non-equal length")
#   if (any(is.na(param) | param > upper | param < lower))
#     return("Parameter value out of bound")
#   else return (TRUE)
# }
###################################################
##### copula class
###################################################
setClass("copula", 
         representation(dimension = "numeric",
                        parameters = "numeric",
                        param.names = "character",
                        param.lowbnd = "numeric",
                        param.upbnd = "numeric",
                        message = "character"),
         #validity = validCopula,
         validity = function(object) {
           dim <- object@dimension
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
         contains = list()
         )


copula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("normal", "t", "clayton", "frank", "gumbel")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  if (fam <= 2) ellipCopula(family, param, dim, ...)
  else archmCopula(family, param, dim, ...)
}


#####################################################
#### show, plot, methods
#####################################################
showCopula <- function(object) {
  cat(object@message, ".\n")
  cat("Dimension: ", object@dimension, "\n")
  cat("Parameters:\n")
  for (i in (1:length(object@parameters)))
    cat("  ", object@param.names[i], " = ", object@parameters[i], "\n")
}

setMethod("show", signature("copula"), showCopula)

#####################################################
####### general methods for all copulas 
#####################################################
dcopula <- function(copula, u) {
  UseMethod("dcopula")
}

pcopula <- function(copula, u) {
  UseMethod("pcopula")
}

rcopula <- function(copula, n) {
  UseMethod("rcopula")
}

kendallsTau <- function(copula, ...) {
  ## bivariate association measurement
  UseMethod("kendallsTau")
}

spearmansRho <- function(copula, ...) {
  ## bivariate association measurement
  UseMethod("spearmansRho")
}

tailIndex <- function(copula, ...) {
  ## bivariate association measurement
  UseMethod("tailIndex")
}

calibKendallsTau <- function(copula, tau) {
  UseMethod("calibKendallsTau")
}

calibSpearmansRho <- function(copula, rho) {
  UseMethod("calibSpearmansRho")
}

#### numerical computation of association measures
kendallsTauCopula <- function(copula, ...) {
  integrand <- function(u) pcopula(copula, u) * dcopula(copula, u)
  .eps <- .Machine$double.eps^0.9
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  4 * integ - 1
}

spearmansRhoCopula <- function(copula, ...) {
  integrand <- function(u) pcopula(copula, u)
  .eps <- .Machine$double.eps^0.9
  lower <- c(.eps, .eps)
  upper <- c(1 - .eps, 1 - .eps)
  integ <- adapt(ndim=2, lower=lower, upper=upper, functn=integrand, ...)$value
  12 * integ - 3                 
}


#### numerical tail index, not accurate
tailIndexCopula <- function(copula, eps = .Machine$double.eps^0.5) {
  u <- eps
  v <- 1 - u
  lower <- pcopula(copula, c(u, u))/u
  upper <- (1 - 2 * v + pcopula(copula, c(v, v)))/ (1 - v)
  c(lower=lower, upper=upper)
}

setMethod("kendallsTau", signature("copula"), kendallsTauCopula)
setMethod("spearmansRho", signature("copula"), spearmansRhoCopula)
setMethod("tailIndex", signature("copula"), tailIndexCopula)

#### numerical calibration
calibKendallsTauCopula <- function(copula, tau) {
  myfun <- function(theta) {
    copula@parameters <- theta
    kendallsTau(copula) - tau
  }
  .eps <- .Machine$double.eps^.5
  lower <- copula@param.lowbnd + .eps
  upper <- copula@param.upbnd - .eps
  lower <- ifelse(lower == -Inf, -.Machine$double.xmax^.5, lower)
  upper <- ifelse(upper == Inf, .Machine$double.xmax^.5, upper)
  sol <- uniroot(myfun, interval=c(lower, upper))$root
  sol
}

calibSpearmansRhoCopula <- function(copula, rho) {
  myfun <- function(theta) {
    copula@parameters <- theta
    spearmansRho(copula) - rho
  }
  .eps <- .Machine$double.eps^.5
  lower <- copula@param.lowbnd
  upper <- copula@param.upbnd
  lower <- ifelse(lower == -Inf, -.Machine$double.xmax^.5, lower)
  upper <- ifelse(upper == Inf, .Machine$double.xmax^.5, upper)
  sol <- uniroot(myfun, interval=c(lower, upper))$root
  sol
}

setMethod("calibKendallsTau", signature("copula"), calibKendallsTauCopula)
setMethod("calibSpearmansRho", signature("copula"), calibSpearmansRhoCopula)

###############################################################
#### elliptical copulas, contains normalCopula and tCopula
###############################################################
validRho <- function(dispstr, dim, lenRho) {
  if (dispstr == "ar1" || dispstr == "ex")
    if (lenRho != 1) return ("Param should have length 1 for dispstr == ar1 or ex")
  if (dispstr == "un")
    if (lenRho != dim * (dim - 1) / 2)
      return("Param should have length dim * (dim - 1) / 2 for dispstr == un")
  if (dispstr == "toep")
    if (lenRho != dim - 1)
      return("Param should have length dim - 1 for dispstr == toep")
  return(TRUE)
}

validEllipCopula <- function(object) {
  dispstr <- object@dispstr
  if (is.na(match(dispstr, c("ar1", "ex", "toep", "un"))))
    return ("dispstr not supported")
  dim <- object@dimension
  rho <- object@getRho(object)
  validRho(dispstr, dim, length(rho))
}

setClass("ellipCopula",
         representation = representation("copula",
           dispstr = "character", getRho="function"),
         validity = validEllipCopula,
         contains = list("copula")
         )

getSigma <- function(copula) {
  dim <- copula@dimension
  rho <- copula@getRho(copula)
  sigma <- diag(dim)
  if (copula@dispstr == "ex") {
    sigma[lower.tri(sigma)] <- rho[1]
    sigma[upper.tri(sigma)] <- rho[1]
  }
  else if (copula@dispstr == "ar1") {
    for (i in 1:dim)  for (j in 1:dim)  sigma[i,j] <- rho ^ abs(i - j)
  }
  else if (copula@dispstr == "un") {
    sigma[lower.tri(sigma)] <- rho
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  }
  else if (copula@dispstr == "toep") {
    for (i in 1:dim) for (j in 1:dim)
      if (i != j) sigma[i,j] <- rho[abs(i - j)]
  }
  sigma
}

ellipCopula <- function(family, param, dim = 2, dispstr = "ex", df = 5, ...) {
  familiesImplemented <- c("normal", "t")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   normalCopula(param, dim = dim, dispstr = dispstr),
                   tCopula(param, dim = dim, dispstr = dispstr, df = df)
                   )
  copula
}

calibKendallsTauEllipCopula <- function(copula, tau) {
  sin((tau * pi) / 2)
}

calibSpearmansRhoEllipCopula <- function(copula, rho) {
  sin (pi * rho / 6) * 2
}

setMethod("calibKendallsTau", signature("ellipCopula"), calibKendallsTauEllipCopula)
setMethod("calibSpearmansRho", signature("ellipCopula"), calibSpearmansRhoEllipCopula)

############################################################
##### archimedean copulas, contains gumbel, frank, ...
############################################################
setClass("archmCopula",
         representation = representation("copula",
           exprdist = "expression"),
         contains = list("copula")
         )

archmCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("clayton", "frank", "gumbel", "amh")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   claytonCopula(param, dim = dim),
                   frankCopula(param, dim = dim),
                   gumbelCopula(param, dim = dim),
                   amhCopula(param, dim = 2)
                   )
  copula
}

genFun <- function(copula, u) {
  UseMethod("genFun")
}

genInv <- function(copula, s) {
  UseMethod("genInv")
}


genFunDer1 <- function(copula, u) {
  UseMethod("genFunDer1")
}

genFunDer2 <- function(copula, u) {
  UseMethod("genFunDer2")
}

## genFunDer <- function(copula, u, n) {## nth derivative
##   UseMethod("genFunDer")
## }

## genInvDer <- function(copula, s, n) {## nth derivative
##   UseMethod("genInvDer")
## }

kendallsTauArchmCopula <- function(copula) {
  integrand <- function(x) genFun(copula, x) / genFunDer1(copula, x)
  1 + 4 * integrate(integrand, 0, 1)$value
}

setMethod("kendallsTau", signature("archmCopula"), kendallsTauArchmCopula)

#######################################################
#### extreme value copulas
#######################################################

setClass("evCopula",
         representation = representation("copula"),
         contains = list("copula")
         )

evCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("galambos", "gumbel", "huslerReiss")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   galambosCopula(param),
                   gumbelCopula(param),
                   huslerReissCopula(param)
                   )
  copula
}
  

Afun <- function(copula, w) {
  UseMethod("Afun")
}

tailIndexEvCopula <- function(copula) {
  lower <- 0
  upper <- 2 - 2 * Afun(copula, 0.5)
  c(lower=lower, upper=upper)
}

setMethod("tailIndex", signature("evCopula"), tailIndexEvCopula)

#######################################################
#### multivariate distibution via copula
#######################################################
setClass("mvdc",
         representation(copula = "copula",
                        margins = "character",
                        paramMargins = "list")
         )

mvdc <- function(copula, margins, paramMargins) {
  val <- new("mvdc", copula = copula, margins = margins, paramMargins = paramMargins)
}
