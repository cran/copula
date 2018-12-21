### ../DESCRIPTION -- Collate:  want setClass() and setGeneric() first
### ~~~~~~~~~~~~~~    --------

## as we have ./Classes.R  ./AllClass.R  before this,
## currently, many setGeneric()s are there

##' A  (generic!) Function with methods to replace the 'fullname' slot :
setGeneric("describeCop", function(x, kind = c("short", "very short", "long"),
                                   prefix = "", ...)
    standardGeneric("describeCop"), signature = c("x", "kind"))

## This could go back to ./dC-dc.R  if all "dCdu" methods also go there
##                       ~~~~~~~~~
setGeneric("dCdu", function(copula, u, ...) standardGeneric("dCdu"))
setGeneric("dCdtheta", function(copula, u, ...) standardGeneric("dCdtheta"))
setGeneric("dlogcdu", function(copula, u, ...) standardGeneric("dlogcdu"))
setGeneric("dlogcdtheta", function(copula, u, ...) standardGeneric("dlogcdtheta"))



setGeneric("fitCopula", function(copula, data, ...) standardGeneric("fitCopula"))

setGeneric("gofCopula", function(copula, x, ...) standardGeneric("gofCopula"))

## parameter names of _free_ parameters:
setGeneric("paramNames", function(x) standardGeneric("paramNames"))
setMethod("paramNames", "xcopula", function(x) paramNames(x@copula))

##' @title Setting the parameter in a copula/Copula
##' @param x (general) copula
##' @param value parameter value
##' @param na.ok logical indicating if NA values are ok for theta
##' @param noCheck logical indicating if parameter constraints should be checked
##' @return a copula of the same class as \code{x} with theta set to \code{value}
##' @author Martin Maechler
setGeneric("setTheta", function(x, value,
                                na.ok = TRUE, noCheck = FALSE, freeOnly=TRUE, ...)
    standardGeneric("setTheta"), signature = c("x", "value"))

setMethod("setTheta", "xcopula", ## set parameter for  daughter copula
	  function(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly=TRUE, ...) {
	      x@copula <- setTheta(x@copula, value, na.ok=na.ok, noCheck=noCheck,
                                   freeOnly=freeOnly, ...)
	      x
	  })


## get parameters
setGeneric("getTheta", function(copula, freeOnly = TRUE, attr = FALSE, named = attr)
    standardGeneric("getTheta"), signature = "copula")
## -- general methods --
setMethod("getTheta", "xcopula", ## from  daughter copula
	  function(copula, freeOnly = TRUE, attr = FALSE, named = attr)
	      getTheta(copula@copula, freeOnly=freeOnly, attr=attr, named=named))
##' default methods [e.g. for khoudraj & indepCopula]: empty parameter
setMethod("getTheta", "parCopula", function(copula) numeric())


## assign free parameter *value*s:
## __NOT_exported__ <= "inconsistent" with 'fixedParam<-'
setGeneric("freeParam<-", function(copula, value) standardGeneric("freeParam<-"))
## set or modify "fixedness" of parameters (TRUE/FALSE) __exported__
setGeneric("fixedParam<-", function(copula, value) standardGeneric("fixedParam<-"))
## logical vector indicating which parameters are free
setGeneric("isFree", function(copula) standardGeneric("isFree"))

## number of (all/free) parameters
setGeneric("nParam", function(copula, freeOnly = FALSE) standardGeneric("nParam"),
           signature = "copula")

##' default methods: empty parameter
setMethod("isFree",   "parCopula", function(copula) logical())
setMethod("nParam",   "parCopula", function(copula) 0L)

