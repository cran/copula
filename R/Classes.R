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


## NB: Mother classes "Copula" (and "parCopula") are in ./AllClass.R
##							  ~~~~~~~~~~

## The independence copula; don't want it to inherit "funny slots" =>
setClass("indepCopula", contains = c("dimCopula", "parCopula"),
         slots = c(exprdist = "expression"))
# no validity needed anymore: have no 'parameters'

### Basic copula class #########################################################

setClass("copula", contains = c("dimCopula", "parCopula", "VIRTUAL"),
         slots = c(parameters = "numeric",
                   param.names = "character",
                   param.lowbnd = "numeric",
                   param.upbnd = "numeric",
                   ## TODO: "a vector" of intervals" paraInterval = "maybeInterval", # [.,.]  (.,.], etc ..
                   fullname = "character"), ## <- DEPRECATED for describeCop(): see ../TODO
         prototype = prototype(dimension = 2L, parameters = NA_real_),
         validity = ##' Check validity of "copula"
         function(object) {
	     dim <- object@dimension # "integer" by definition
	     if (length(dim) != 1L) return("'dim' must be an integer (>= 2)")
	     if (dim < 2) return("dim must be >= 2")
             param <- object@parameters
             upper <- object@param.upbnd
             lower <- object@param.lowbnd
             lp <- length(param) # if that is > 1, we allow |upper| or |lower| = 1
             if(!(lp == length(upper) || (lp > 1L && length(upper) == 1L)))
		 return("Parameter and upper bound have non-equal length")
             if(!(lp == length(lower) || (lp > 1L && length(lower) == 1L)))
		 return("Parameter and lower bound have non-equal length")
             intervChk <- ## TODO: mkParaConstr(object@paraInterval)
                 function(par) all(is.na(param) | (lower <= param & param <= upper))
             ina.p <- is.na(param)
             if(!all(ina.p)) {
		 ##if(any(ina.p)) return("some (but not all) parameter values are  NA")
                 if(!intervChk(param)) return("Parameter value(s) out of bound")
             }
	     ## to allow (all) NA parameters: "not-yet-parametrized"
	     TRUE
         })


## Methods for 'copula'
setGeneric("dCopula", function(u, copula, log=FALSE, ...) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    stopifnot(dim(copula) == ncol(u))
    u.is.out <- outside.01(u, strictly=FALSE)## on.boundary _or_ outside
    if(any.out <- any(u.is.out, na.rm=TRUE))
	u[] <- pmax(0, pmin(1, u)) # <- "needed", as some methods give error
    r <- standardGeneric("dCopula") # the result of calling  <dCopula-method)(u, copula, ..)
    if(any.out) ## on boundary _or_ outside cube  ==> zero mass :
	r[u.is.out & !is.na(u.is.out)] <- if(log) -Inf else 0.
    r
})
setGeneric("pCopula", function(u, copula, ...) {
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    ## here as well, 'outside' and 'on-boundary' are equivalent:
    u[] <- pmax(0, pmin(1, u))
    standardGeneric("pCopula")
})
setGeneric("rCopula", function(n, copula, ...) standardGeneric("rCopula"))
setGeneric("tau", function(copula, ...) standardGeneric("tau"))
setGeneric("rho", function(copula, ...) standardGeneric("rho"))
setGeneric("lambda", function(copula, ...) standardGeneric("lambda"))
setGeneric("iTau", function(copula, tau, ...) standardGeneric("iTau"))
setGeneric("iRho", function(copula, rho, ...) standardGeneric("iRho"))
setGeneric("dTau", function(copula, ...) standardGeneric("dTau"))
setGeneric("dRho", function(copula, ...) standardGeneric("dRho"))
setGeneric("dTauFun", function(copula) standardGeneric("dTauFun"))
setGeneric("dRhoFun", function(copula) standardGeneric("dRhoFun"))

## Defunct methods
calibKendallsTau <- function(copula, tau) { .Defunct("iTau"); iTau(copula,tau) }
calibSpearmansRho <- function(copula, rho) { .Defunct("iRho"); iRho(copula,rho) }
kendallsTau <- function(copula) { .Defunct("tau"); tau(copula) }
spearmansRho <- function(copula) { .Defunct("rho"); rho(copula) }
genInv <- function(copula, s) { .Defunct("psi"); psi(copula,s) }
genFun <- function(copula, u) { .Defunct("iPsi"); iPsi(copula, u) }
genFunDer1 <- function(copula, u){ .Defunct("diPsi"); diPsi(copula, u) }
genFunDer2 <- function(copula, u){ .Defunct("diPsi(*, degree=2)"); diPsi(copula, u, degree=2) }
AfunDer <- function(copula, w) { .Defunct("dAdu"); dAdu(copula, w) }
Afun    <- function(copula, w) { .Defunct("A"); A(copula, w) }
pcopula <- function(copula, u, ...) { .Defunct("pCopula"); pCopula(u, copula) }
dcopula <- function(copula, u, ...) { .Defunct("dCopula"); dCopula(u, copula, ...) }
rcopula <- function(copula, n, ...) { .Defunct("rCopula"); rCopula(n, copula, ...) }

## Deprecated methods
tailIndex <- function(copula) { .Deprecated("lambda"); lambda(copula) } # mid 2016


### Frechet--Hoeffding bounds ##################################################

setClass("fhCopula", contains = c("dimCopula", "parCopula", "VIRTUAL"),
         slots = c(exprdist = "expression")) # -> used in dC-dc.R

## Lower Frechet--Hoeffding bound W
setClass("lowfhCopula", contains = "fhCopula")

## Upper Frechet--Hoeffding bound M
setClass("upfhCopula", contains = "fhCopula")

##-> fhCopula.R, lowfhCopula.R, upfhCopula.R  for implementation
##   ~~~~~~~~~~  ~~~~~~~~~~~~~  ~~~~~~~~~~~~

### Empirical copula ###########################################################

setClass("empCopula", contains = "Copula",
         slots = c(X = "matrix",
		   smoothing = "character",
		   offset = "numeric",
		   ties.method = "character"),
	 ## validity methods must return TRUE or string
	 validity = function(object) {
	    X <- object@X
	    if(ncol(X) < 2L)
		"X must have at least 2 columns"
	    else if(anyNA(X) || any(X < 0) || any(X > 1))
		"All entries of X[] must be in [0,1]"
	    else
		TRUE
	 })


### Elliptical copulas (normalCopula, tCopula) #################################

## NOTE: This is related to npar.ellip() in ./ellipCopula.R
validRho <- function(dispstr, dim, lenRho) {
    switch(dispstr, ## checking for correct 'dispstr'
	   "ar1" =, "ex" = {
	       if (lenRho != 1)
		   return(gettextf("'rho' parameter should have length 1 for 'dispstr' = \"%s\"",
				   dispstr))
	   },
	   "un" = {
	       if (lenRho != (L <- dim * (dim - 1) / 2))
		   return(gettextf("'rho' parameter should have length dim * (dim - 1) / 2 (= %d) for 'dispstr' = \"%s\"",
				   L, dispstr))
	   },
	   "toep" = {
	       if (lenRho != dim - 1)
		   return(gettextf("'rho' parameter should have length dim-1 (= %d) for 'dispstr' = \"%s\"",
				   dim-1, dispstr))
	   },
	   ## otherwise
	   return("'dispstr' not supported (yet)"))

    TRUE
}

setClass("ellipCopula", contains = c("copula", "VIRTUAL"),
	 slots = c(dispstr = "character", getRho = "function"),
         validity = function(object)
             validRho(dispstr=object@dispstr, dim=object@dimension,
                      length(object@getRho(object))))

setClass("normalCopula", contains = "ellipCopula"
         ## not really needed -- validity for ellipCopula is checked already
         ## , validity =  function(object) {
         ##     can do more if needed here
         ## }
         )

setClass("tCopula", contains = "ellipCopula",
	 , slots = c(df.fixed = "logical")
         ## validity =  function(object) {
         ##     ## JY: making sure @df == tail(.@parameter, 1)
         ##     if (has.par.df(object))
         ##         stopifnot(object@df == tail(object@parameters, 1))
         ##     TRUE
         ## }
         )


### Archimedean copulas (AMH, Clayton, Frank, Gumbel, Joe) #####################

setClass("archmCopula", contains = c("copula", "VIRTUAL"),
         slots = c(exprdist = "expression"),
         ## extra check in addition to "copula" validity:
         validity = function(object) {
	     if(length(object@parameters) != 1)
		 "Our Archimedean copulas must currently have a 1-dimensional parameter"
	     else
		 TRUE
         })

## Clayton
setClass("claytonCopula", contains = "archmCopula")

## Gumbel (also an evCopula) --> see below

## Frank
setClass("frankCopula", contains = "archmCopula")

## AMH
setClass("amhCopula", contains = "archmCopula")

## Joe
setClass("joeCopula", contains = "archmCopula")

## These can have negative tau for d = 2 only:
## setClassUnion(..) # <- over kill?; should be enough:
archm.neg.tau <- c("amhCopula", "claytonCopula", "frankCopula")

## Methods for 'archmCopula'
setGeneric("psi", function(copula, s) standardGeneric("psi"))
##FIXME 'log' compulsory:
##setGeneric("iPsi", function(copula, u, log, ...) standardGeneric("iPsi"))
setGeneric("iPsi", function(copula, u, ...) standardGeneric("iPsi"))
setGeneric("dPsi", function(copula, s, ...) standardGeneric("dPsi"))
setGeneric("diPsi", function(copula, u, degree=1, log=FALSE, ...) standardGeneric("diPsi"))


### Extreme value copulas (Galambos, Husler--Reiss, Gumbel, Tawn, tEV) #########

setClass("evCopula", contains = c("copula", "VIRTUAL"))

## Galambos
setClass("galambosCopula", contains = "evCopula",
	 slots = c(exprdist = "expression"))

## Gumbel (also an Archimedean copula)
setClass("gumbelCopula", contains = list("archmCopula", "evCopula"))

## Husler--Reiss
setClass("huslerReissCopula", contains = "evCopula",
	 slots = c(exprdist = "expression"))

## Tawn (does not offer full range of dependence)
setClass("tawnCopula", contains = "evCopula",
	 slots = c(exprdist = "expression"))

## tEV
setClass("tevCopula", contains = "evCopula"
	 , slots = c(df.fixed = "logical")
         ## , validity = function(object) {
         ##     ## JY: making sure @df == tail(@parametersnot, 1)
         ##     if(object@df != tail(object@parameters, 1))
         ##         return("'df' must be the last element of @parameters")
         ##     TRUE
         ## }
         )

setGeneric("A", function(copula, w) standardGeneric("A"))
setGeneric("dAdu", function(copula, w) standardGeneric("dAdu"))
setGeneric("dAdtheta", function(copula, w) standardGeneric("dAdtheta"))


### Other copulas ##############################################################

## Marshall--Olkin copula
setClass("moCopula", contains = "copula", slots = c(exprdist = "expression"))

## Farlie-Gumbel-Morgenstern copula
setClass("fgmCopula", contains = "copula", slots = c(exprdist = "expression",
                                                     subsets.char = "character"),
         ## verify that the pdf is positive at each vertex of [0,1]^dim
	 validity = function(object) {
	     dim <- object@dimension
	     if (dim < 2)
		 "must be 2-dimensional"
	     else if (dim == 2)
		 TRUE
	     else { # dim >= 3
		 param <- object@parameters
		 valid <- .C(validity_fgm,
			     as.integer(dim),
			     as.double(c(rep(0,dim+1),param)),
			     valid = integer(1))$valid
		 if (valid == 0) "Bad vector of parameters" else TRUE
	     }
	 })

## Plackett copula
setClass("plackettCopula", contains = "copula", slots = c(exprdist = "expression"))


### Multivariate distibutions via copulas ######################################

validMvdc <- function(object) {
    dim <- dim(object@copula) # @dimension may not exist
    if(!is.finite(dim) || dim < 2)
	return("'dimension' must be integer >= 2")
    if(dim != length(object@margins))
	return("'dimension' does not match margins' length")
    if(dim != length(pm <- object@paramMargins))
	return("'dimension' does not match paraMargins' length")
    if(!all(vapply(pm, function(e) is.list(e) || is.numeric(e), NA)))
	return("'paramMargins' elements must all be list()s or numeric vectors")
    okNms <- function(nms) !is.null(nms) && all(nzchar(nms))
    if(object@marginsIdentical) {
	if(!all(object@margins[1] == object@margins[-1]))
	    return("margins are not identical")
	pm1 <- pm[[1]]
	for(i in 2:dim) {
	    if(!identical(pm1, pm[[i]]))
		return("margins are not identical")
	}
	if(length(pm1) > 0 && !okNms(names(pm1)))
	    return("'paramMargins' must be named properly")
    }
    else ## not identical margins: check each
	for(i in seq_len(dim)) {
	    pmi <- pm[[i]]
	    if(length(pmi) > 0 && !okNms(names(pmi)))
		return(gettextf("'paramMargins[[%d]]' must be named properly", i))
	    ## TODO(?): check more similar to (/ instead of) those in mvdc() --> ./mvdc.R
	}
    TRUE
}## validMvdc()

setClass("mvdc", contains = "xcopula", #-> slot "copula"; additionally :
	 slots = c(
                   margins = "character",
                   paramMargins = "list",
                   marginsIdentical = "logical"),
	 validity = validMvdc)

## methods like {dpr}mvdc are defined in mvdc.R

## A fitted multivariate distribution -- "generic mother class",
## "fitCopula" and "fitMvdc" will inherit from it:
setClass("fittedMV",
	 slots = c(estimate = "numeric",
		   var.est = "matrix", ##  and matrix(,0,0) means "not estimated/available"
                   loglik = "numeric",
                   nsample = "integer",
                   method = "character",
                   call = "call",
                   ## convergence = "integer",
                   fitting.stats = "list"))

coef.fittedMV <- function(object, SE = FALSE, ...) {
    pNms <- paramNames(object)
    if(SE)
	structure(cbind(object@estimate,
			if(length(V <- object@var.est)) sqrt(diag(V))
			else rep(NA_real_, length(pNms)),
			deparse.level=0L),
		  dimnames = list(pNms, c("Estimate", "Std. Error")))
    else
	setNames(object@estimate, pNms)
}

nobs.fittedMV <- function(object, ...) object@nsample

vcov.fittedMV <- function(object, ...) {
    pNms <- paramNames(object)
    structure(if(length(V <- object@var.est)) V
              else matrix(NA, length(pNms), length(pNms)),
              dimnames = list(pNms, pNms))
}

logLik.fittedMV <- function(object, ...) {
    structure(object@loglik,
	      nobs = object@nsample,
	      df = length(object@estimate), # == #{free param.}
	      class = "logLik")
}

##' Does the copula have 'df' as *free*, i.e., non-fixed parameter?
has.par.df <- function(cop, classDef = getClass(class(cop)),
                       isEllip = extends(classDef, "ellipCopula")) {
    ((isEllip && extends(classDef, "tCopula")) || extends(classDef, "tevCopula")) &&
        !cop@df.fixed
}
