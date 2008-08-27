###################################################
##### basic copula class
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

## general methods for copula
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
  UseMethod("kendallsTau")
}

spearmansRho <- function(copula, ...) {
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

tauDer <- function(copula) {
  UseMethod("tauDer")
}

rhoDer <- function(copula) {
  UseMethod("rhoDer")
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

## normal copula
validNormalCopula <- function(object) {
  validEllipCopula(object)
  ## can do more if needed here
}

setClass("normalCopula",
         representation = representation("ellipCopula"),
         validity = validNormalCopula,
         contains = list("copula", "ellipCopula")
         )

## t copula
validTCopula <- function(object) {
  df <- getdf(object)
  if (df <= 0) return ("df should be > 0")
  validEllipCopula(object)
}

setClass("tCopula",
         representation = representation("ellipCopula",
           df = "numeric",
           df.fixed = "logical"),
         validity = validTCopula,
         contains = list("copula", "ellipCopula")
         )


## methods for ellipCopula??


############################################################
#### archimedean copulas, contains clayton, gumbel, frank,
#### amh, ...
############################################################
setClass("archmCopula",
         representation = representation("copula",
           exprdist = "expression"),
         contains = list("copula")
         )

## clayton copula
setClass("claytonCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )

## gumbel copula, also an ev copula

## frank copula
setClass("frankCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )

## amh copula
setClass("amhCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula")
         )

## methods for archmCopulas
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

#######################################################
#### extreme value copulas, contains galambos, husler-reiss,
#### gumbel, ...
#######################################################

setClass("evCopula",
         representation = representation("copula"),
         contains = list("copula")
         )

## galambos copula
setClass("galambosCopula",
         representation = representation("evCopula"),
         contains = list("copula", "evCopula")
         )

## gumbel copula, also an archm copula; how to clean this up?
setClass("gumbelCopula",
         representation = representation("archmCopula"),
         contains = list("copula", "archmCopula", "evCopula")
         )

## husler-reiss copula
setClass("huslerReissCopula",
         representation = representation("evCopula"),
         contains = list("copula", "evCopula")
         )


## methods for evCopula
Afun <- function(copula, w) {
  UseMethod("Afun")
}


#######################################################
#### other copulas
#######################################################

## Farlie-Gumbel-Morgenstern multivariate copula
setClass("fgmCopula",
         representation = representation("copula",
         exprdist = "expression"),
         ## verify that the pdf is positive at each vertex of [0,1]^dim
         validity = function(object) {
             dim <- object@dimension
             if (dim == 2)
                 return(TRUE);
             param <- object@parameters
             valid <- .C("validity_fgm", 
                         as.integer(dim),
                         as.double(c(rep(0,dim+1),param)),
                         valid = integer(1),
                         PACKAGE="copula")$valid
             if (valid == 0)
                 return("Bad vector of parameters")
             else
                 return(TRUE)   
         },
         contains = list("copula")
         )

## plackett copula
setClass("plackettCopula",
         representation = representation("copula",
           exprdist = "expression"),
         contains = list("copula")
         )

#######################################################
#### multivariate distibution via copula
#######################################################
setClass("mvdc",
         representation(copula = "copula",
                        margins = "character",
                        paramMargins = "list")
         )

## methods like {dpr}mvdc are defined in mvdc.R
