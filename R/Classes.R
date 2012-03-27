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


### basic copula class #########################################################

##' Check validity of "copula"  (not exported for now)
validCopula <- function(object) {
    dim <- object@dimension
    if (dim != as.integer(dim))
        return("dim must be integer")
    if (dim < 2)
        return("dim must be >= 2")
    param <- object@parameters
    upper <- object@param.upbnd
    lower <- object@param.lowbnd
    lp <- length(param)
    if (lp != length(upper) && length(upper) != 1)
        return("Parameter and upper bound have non-equal length")
    if (lp != length(lower) && length(lower) != 1)
        return("Parameter and lower bound have non-equal length")
    if (any(is.na(param) | param > upper | param < lower))
        return("Parameter value out of bound")
    else return (TRUE)
}

setClass("copula",
         representation(dimension = "integer", # as for "nacopula"
                        parameters = "numeric",
                        param.names = "character",
                        param.lowbnd = "numeric",
                        param.upbnd = "numeric",
                        message = "character",
                        "VIRTUAL"),
         validity = validCopula)

## general methods for copula
setGeneric("dcopula", function(copula, u, log=FALSE, ...) standardGeneric("dcopula"))
setGeneric("pcopula", function(copula, u, ...) standardGeneric("pcopula"))
setGeneric("rcopula", function(copula, n, ...) standardGeneric("rcopula"))
setGeneric("kendallsTau", function(copula) standardGeneric("kendallsTau"))
setGeneric("spearmansRho", function(copula) standardGeneric("spearmansRho"))
setGeneric("tailIndex", function(copula, ...) standardGeneric("tailIndex"))
setGeneric("calibKendallsTau", function(copula, tau) standardGeneric("calibKendallsTau"))
setGeneric("calibSpearmansRho", function(copula, rho) standardGeneric("calibSpearmansRho"))
setGeneric("tauDer", function(copula, ...) standardGeneric("tauDer"))
setGeneric("rhoDer", function(copula, ...) standardGeneric("rhoDer"))

setGeneric("tauDerFun", function(copula) standardGeneric("tauDerFun"))
setGeneric("rhoDerFun", function(copula) standardGeneric("rhoDerFun"))


### independent copula class ###################################################

setClass("indepCopula", contains = "copula",
         representation(exprdist = "expression"))

### elliptical copulas, contains normalCopula and tCopula ######################

validRho <- function(dispstr, dim, lenRho) {
    switch(dispstr, ## checking for correct 'dispstr'
	   "ar1" =, "ex" = {
	       if (lenRho != 1)
		   return(sprintf("'rho' parameter should have length 1 for 'dispstr' = \"%s\"",
				  dispstr))
	   },
	   "un" = {
	       if (lenRho != dim * (dim - 1) / 2)
		   return("Param should have length dim * (dim - 1) / 2 for dispstr == un")
	   },
	   "toep" = {
	       if (lenRho != dim - 1)
		   return("Param should have length dim - 1 for dispstr == toep")
	   },
	   ## otherwise
	   return("'dispstr' not supported (yet)"))

    return(TRUE)
}

validEllipCopula <- function(object) {
  rho <- object@getRho(object)
  validRho(dispstr=object@dispstr, dim=object@dimension,
           length(rho))
}

setClass("ellipCopula", contains = "copula",
         representation(dispstr = "character", getRho="function"),
         validity = validEllipCopula)

if(FALSE) # not yet needed -- validEllipCopula() is used anyway
##' normal copula
validNormalCopula <- function(object) {
  validEllipCopula(object)
  ## can do more if needed here
}
setClass("normalCopula", contains = "ellipCopula")
         ## validity = validNormalCopula)

## t copula
validTCopula <- function(object) {
  df <- getdf(object)
  if (df <= 0) return ("df should be > 0")
  validEllipCopula(object)
}

setClass("tCopula", representation(df = "numeric", df.fixed = "logical"),
         contains = "ellipCopula",
         validity = validTCopula)


## methods for ellipCopula??


### Archimedean copulas, contains AMH, Clayton, Frank, Gumbel, ... #############

setClass("archmCopula", representation(exprdist = "expression", "VIRTUAL"),
	 contains = "copula")

## clayton copula
setClass("claytonCopula", contains = "archmCopula")

## gumbel copula, also an ev copula

## frank copula
setClass("frankCopula", contains = "archmCopula")

## amh copula
setClass("amhCopula", contains = "archmCopula")

## methods for archmCopulas
setGeneric("genFun", function(copula, u) standardGeneric("genFun"))
setGeneric("genInv", function(copula, s) standardGeneric("genInv"))
setGeneric("genFunDer1", function(copula, u) standardGeneric("genFunDer1"))
setGeneric("genFunDer2", function(copula, u) standardGeneric("genFunDer2"))


### Extreme value copulas, contains galambos, husler-reiss, gumbel, ... ########

setClass("evCopula", representation("VIRTUAL"), contains = "copula")

## galambos copula
setClass("galambosCopula", representation(exprdist = "expression"),
         contains = "evCopula")

## gumbel copula, also an archm copula;
setClass("gumbelCopula", contains = list("archmCopula", "evCopula"))

## husler-reiss copula
setClass("huslerReissCopula",representation(exprdist = "expression"),
         contains = "evCopula")

## tawn copula; does not offer full range of dependence
setClass("tawnCopula", representation(exprdist = "expression"),
         contains = "evCopula")

## tEV copula
setClass("tevCopula", representation(df = "numeric", df.fixed = "logical"),
         contains = "evCopula")

setGeneric("Afun", function(copula, w) standardGeneric("Afun"))
setGeneric("AfunDer", function(copula, w) standardGeneric("AfunDer"))
setGeneric("derAfunWrtParam", function(copula, w) standardGeneric("derAfunWrtParam"))


### Other copulas ##############################################################

## Farlie-Gumbel-Morgenstern multivariate copula
setClass("fgmCopula", representation(exprdist = "expression"),
         contains = "copula",
         ## verify that the pdf is positive at each vertex of [0,1]^dim
         validity = function(object) {
             dim <- object@dimension
             if (dim == 2)
                 return(TRUE)
             param <- object@parameters
             valid <- .C(validity_fgm,
                         as.integer(dim),
                         as.double(c(rep(0,dim+1),param)),
                         valid = integer(1))$valid
             if (valid == 0)
                 return("Bad vector of parameters")
             else
                 return(TRUE)
         })


## plackett copula
setClass("plackettCopula",representation(exprdist = "expression"),
         contains = "copula")

### Multivariate distibution via copula ########################################

setClass("mvdc",
         representation(copula = "copula",
                        margins = "character",
                        paramMargins = "list",
         		marginsIdentical = "logical"),
         validity = function(object) {
           dim <- object@copula@dimension # guaranteed to be >= 2
           if(dim != length(object@margins))
               return("'dimension' does not match margins' length")
           if(dim != length(object@paramMargins))
               return("'dimension' does not match paraMargins' length")
           if(object@marginsIdentical){
             if(!all(object@margins[1] == object@margins[-1]))
               return("margins are not identical")
             pm1 <- object@paramMargins[[1]]
             for(i in 2:dim) {
               if(!identical(pm1, object@paramMargins[[i]]))
                 return("margins are not identical")
             }
           }
           TRUE
         })

## methods like {dpr}mvdc are defined in mvdc.R

###-------------------------- Glue   "copula" <-> "nacopula"

##' The mother of all copula classes:
setClassUnion("Copula",
              members = c("copula", "nacopula"))
## NB: "acopula" *not* : It has no dimension, is rather a family object
