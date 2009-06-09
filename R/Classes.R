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
setGeneric("dcopula", function(copula, u) standardGeneric("dcopula"))
setGeneric("pcopula", function(copula, u) standardGeneric("pcopula"))
setGeneric("rcopula", function(copula, n) standardGeneric("rcopula"))
setGeneric("kendallsTau", function(copula, ...) standardGeneric("kendallsTau"))
setGeneric("spearmansRho", function(copula, ...) standardGeneric("spearmansRho"))
setGeneric("tailIndex", function(copula, ...) standardGeneric("tailIndex"))
setGeneric("calibKendallsTau", function(copula, tau) standardGeneric("calibKendallsTau"))
setGeneric("calibSpearmansRho", function(copula, rho) standardGeneric("calibSpearmansRho"))
setGeneric("tauDer", function(copula, ...) standardGeneric("tauDer"))
setGeneric("rhoDer", function(copula, ...) standardGeneric("rhoDer"))

setGeneric("tauDerFun", function(copula) standardGeneric("tauDerFun"))
setGeneric("rhoDerFun", function(copula) standardGeneric("rhoDerFun"))

###################################################
#### independent copula class
###################################################
## setClass("indepCopula",
##          representation(dimension = "numeric",
##                         message = "character"),
##          contains = list()
##          )
setClass("indepCopula",
         representation("copula"),
         contains = list("copula")
         )


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
         contains = list("ellipCopula")
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
         contains = list("ellipCopula")
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
         contains = list("archmCopula")
         )

## gumbel copula, also an ev copula

## frank copula
setClass("frankCopula",
         representation = representation("archmCopula"),
         contains = list("archmCopula")
         )

## amh copula
setClass("amhCopula",
         representation = representation("archmCopula"),
         contains = list("archmCopula")
         )

## methods for archmCopulas

setGeneric("genFun", function(copula, u) standardGeneric("genFun"))
setGeneric("genInv", function(copula, s) standardGeneric("genInv"))
setGeneric("genFunDer1", function(copula, u) standardGeneric("genFunDer1"))
setGeneric("genFunDer2", function(copula, u) standardGeneric("genFunDer2"))

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
         representation = representation("evCopula",
           exprdist = "expression"),
         contains = list("evCopula")
         )

## gumbel copula, also an archm copula; 
setClass("gumbelCopula",
         representation = representation("archmCopula"),
         contains = list("archmCopula", "evCopula")
         )

## husler-reiss copula
setClass("huslerReissCopula",
         representation = representation("evCopula",
           exprdist = "expression"),
         contains = list("evCopula")
         )

## tawn copula; does not offer full range of dependence
setClass("tawnCopula",
         representation = representation("evCopula",
           exprdist = "expression"),
         contains = list("evCopula")
         )

## tEV copula
setClass("tevCopula",
         representation = representation("evCopula",
           df = "numeric",
           df.fixed = "logical"),
         contains = list("evCopula")
         )

setGeneric("Afun", function(copula, w) standardGeneric("Afun"))
setGeneric("AfunDer", function(copula, w) standardGeneric("AfunDer"))
setGeneric("derAfunWrtParam", function(copula, w) standardGeneric("derAfunWrtParam"))
           
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
                        paramMargins = "list",
         		marginsIdentical = "logical"),
         validity = function(object){
           dim <- object@copula@dimension
           if(object@marginsIdentical){
             if(!all(object@margins[1] == object@margins[-1]))
               return("margins are not identical")
             for(i in 2:dim){
               if(!identical( object@paramMargins[[1]], object@paramMargins[[i]]))
                 return("margins are not identical")
             }
           }
           TRUE
         }
         )

## methods like {dpr}mvdc are defined in mvdc.R
