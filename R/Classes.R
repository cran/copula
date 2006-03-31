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
####### new general methods for all copulas 
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
                   tCopula(param, dim = dim, dispstr = dispstr, df = df))
  copula
}

##### archimedean copulas, contains gumbel, frank, ...
setClass("archmCopula",
         representation = representation("copula",
           exprdist = "expression"),
         contains = list("copula")
         )

archmCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("clayton", "frank", "gumbel")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   claytonCopula(param, dim = dim),
                   frankCopula(param, dim = dim),
                   gumbelCopula(param, dim = dim))
  copula
}

genFun <- function(copula, u) {
  UseMethod("genFun")
}

genInv <- function(copula, s) {
  UseMethod("genInv")
}


genDer1 <- function(copula, u) {
  UseMethod("genDer1")
}

genDer2 <- function(copula, u) {
  UseMethod("genDer2")
}
##### copula constructor
copula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("normal", "t", "clayton", "frank", "gumbel")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  if (fam <= 2) ellipCopula(family, param, dim, ...)
  else archmCopula(family, param, dim, ...)
}

# #### extreme value copulas
# # setClass("EVCopula",
# #          representation = representation("copula"),
# #          contains = list("copula")
# #          )

# # Afun <- function(copula, u) {
# #   UseMethod("Afun")
# # }

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
