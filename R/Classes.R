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
validEllipCopula <- function(object) {
  dim <- object@dimension
  param <- object@parameters
  ##val <- validCopula(object)
  if (is.na(match(object@corstr, c("ar1", "ex", "toep", "un"))))
    return ("corstr not supported")
  if (object@corstr == "ar1" || object@corstr == "ex")
    if (length(param) != 1) return ("Param should have length 1 for corstr == ar1 or ex")
  if (object@corstr == "un")
    if (length(param) != dim * (dim - 1) / 2)
      return("Param should have length dim * (dim - 1) / 2 for corstr == un")
  if (object@corstr == "toep")
    if (length(param) != dim - 1)
      return("Param should have length dim - 1 for corstr == toep")
  return(TRUE)
}

setClass("ellipCopula",
         representation = representation("copula",
           corstr = "character"),
         validity = validEllipCopula,
         contains = list("copula")
         )

getSigma <- function(copula) {
  dim <- copula@dimension
  param <- copula@parameters
  sigma <- diag(dim)
  if (copula@corstr == "ex") {
    sigma[lower.tri(sigma)] <- param[1]
    sigma[upper.tri(sigma)] <- param[1]
  }
  else if (copula@corstr == "ar1") {
    for (i in 1:dim)  for (j in 1:dim)  sigma[i,j] <- param ^ abs(i - j)
  }
  else if (copula@corstr == "un") {
    sigma[lower.tri(sigma)] <- param
    sigma[upper.tri(sigma)] <- t(sigma)[upper.tri(sigma)]
  }
  else if (copula@corstr == "toep") {
    for (i in 1:dim) for (j in 1:dim)
      if (i != j) sigma[i,j] <- param[abs(i - j)]
  }
  sigma
}

ellipCopula <- function(family, param, dim = 2, corstr = "ex", df = 5, ...) {
  familiesImplemented <- c("normal", "t")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   normalCopula(param, dim = dim, corstr = corstr),
                   tCopula(param, dim = dim, corstr = corstr, df = df))
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
