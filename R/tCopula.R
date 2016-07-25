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


## DONE: allow "param = NA" --
##  2) npar.ellip(dim, dispstr, df.fixed) |-->  "dimension" (length) of param
##  --> TODO: use it ---> define a generic  nparam(.) {with methods for all copula}
### TODO:  {also for "normalCopula"}
##  3) validity should check  pos.definiteness for "un"structured (maybe "toeplitz"
##

##' @title Constructor of tCopula
##' @param param numeric vector of dispersion parameters (df excluded)
##' @param dim integer dimension
##' @param dispstr dispersion structure ("ex", "ar1", "toep", or "un")
##' @param df numeric, degrees of freedom
##' @param df.fixed logical, TRUE = df fixed
##' @return an object of tCopula
tCopula <- function(param = NA_real_, dim = 2L, dispstr = "ex",
		    df = 4, df.fixed = FALSE)
{
    dim <- as.integer(dim)
    if(!is.numeric(param)) storage.mode(param) <- "double" # for NA, keeping attributes!
    stopifnot((pdim <- length(param)) >= 1)
    if(pdim == 1 && is.na(param)) ## extend it (rho only!)
	pdim <- length(param <- rep(param, length.out = npar.ellip(dim, dispstr)))
    parameters <- c(param, df) ## df is another parameter __at end__
    param.names <- c(paste("rho", seq_len(pdim), sep="."), "df")
    param.lowbnd <- c(rep.int(-1, pdim), 1e-6)
    param.upbnd	 <- c(rep.int( 1, pdim), Inf)
    attr(parameters, "fixed") <- c(isFixedP(param), df.fixed)

    new("tCopula",
	dispstr = dispstr,
	dimension = dim,
	parameters = parameters,
	## df = df, ## JY: need to deprecate this in the future
	df.fixed = df.fixed,
	param.names = param.names,
	param.lowbnd = param.lowbnd,
	param.upbnd = param.upbnd,
	fullname = paste("t copula family", if(df.fixed) paste("df fixed at", df)),
	getRho = function(obj) {
	    par <- obj@parameters
            par[-length(par)]
	}
	)
}

### Used for tCopula *and* tevCopula
##' @title Coerce to a Copula with fixed df
##' @param copula a copula object with a "df.fixed" slot
##' @param classDef class definition from getClass
##' @return a copula with df.fixed = TRUE
##' @author MH/MM?
as.df.fixed <- function(copula, classDef = getClass(class(copula))) {
    if (copula@df.fixed) return(copula)
    stopifnot(extends(classDef, "tCopula") || extends(classDef, "tevCopula"))
    if (!is.null(fixed <- attr(copula@parameters, "fixed"))) {
        ## for tCopula with "fixed" attr in @parameters
        fixed[length(fixed)] <- TRUE
    } else {
        fixed <- c(rep(FALSE, length(copula@parameters) - 1L), TRUE)
    }
    attr(copula@parameters, "fixed") <- fixed
    copula@df.fixed <- TRUE # (to be deprecated)
    copula
}

##' @title Get the df parameter of a t/tev copula
##' @param object a copula object
##' @return the df of copula
##' @author Jun Yan
getdf <- function(object) object@parameters[[length(object@parameters)]]

## __ NOWHERE USED __
##
## ##' @title Set the df parameter of a t/tev copula
## ##' @param object a copula object
## ##' @param df the df value to be set
## ##' @return a copula object with df set
## ##' @author Jun Yan
## `setdf<-` <- function(object, value) {
##     stopifnot(is.numeric(value), length(value) == 1, value > 0, value < Inf)
##     object@df <- value
##     if (!is.null(fixed <- attr(object@parameters, "fixed"))) {
##         object@parameters[length(object@parameters)] <- value
##     } else { ## old version
##         fixed <- c(rep(FALSE, length(object@parameters)), object@df.fixed)
##         if (object@df.fixed) {
##             object@parameters <- c(object@parameters, value)
##             object@param.names  <- c(object@param.names, "df")
##             object@param.lowbnd <- c(object@param.lowbnd, 1e-6)
##             object@param.upbnd  <- c(object@param.upbnd, Inf)
##         } else object@parameters[length(object@parameters)] <- value
##         attr(object@parameters, "fixed") <- fixed
##     }
##     object
## }

rtCopula <- function(n, copula) {
  df <- getdf(copula)
  pt(rmvt(n, sigma = getSigma(copula), df = df), df = df)
}


ptCopula <- function(u, copula, ...) {
  dim <- copula@dimension
  i.lower <- rep.int(-Inf, dim)
  sigma <- getSigma(copula)
  df <- getdf(copula)
  if(!(df==Inf || df == as.integer(df)))
    stop("'df' is not integer (or Inf); therefore, pCopula() cannot be computed yet")
  ## more checks now  pCopula() *generic*
  apply(u, 1, function(x) if(any(is.na(x))) NA_real_ else
	pmvt(lower = i.lower, upper = qt(x, df = df), sigma = sigma, df = df, ...))
}

dtCopula <- function(u, copula, log = FALSE, ...) {
  dim <- copula@dimension
  sigma <- getSigma(copula)
  df <- getdf(copula)
  ## more checks now  dCopula() *generic*
  r <- numeric(nrow(u)) # i.e. 0  by default (i.e. "outside")
  ok <- u.in.01(u)
  x <- qt(u[ok, , drop=FALSE], df)
  ## work in log-scale [less over-/under-flow, then (maybe) transform]:
  r[ok] <- dmvt(x, delta = rep.int(0, dim), sigma = sigma, df = df, log = TRUE) -
      rowSums(dt(x, df = df, log=TRUE))
  if(log) r else exp(r)
}

showTCopula <- function(object) {
  print.copula(object)
  if (object@dimension > 2) cat("dispstr: ", object@dispstr, "\n")
  if (object@df.fixed) cat("df is fixed at", getdf(object), "\n")
  invisible(object)
}

lambdaTCopula <- function(copula)
{
  ## McNeil, Frey, Embrechts (p. 211, 2005)
    df <- getdf(copula)
    rho <- copula@getRho(copula)
    if(is.infinite(df)) {
        lambdaNormalCopula(normalCopula(rho))
    } else {
        res <- 2 * pt(- sqrt((df + 1) * (1 - rho) / (1 + rho)), df=df + 1)
        c(lower = res, upper = res)
    }
}

setMethod("rCopula", signature("numeric", "tCopula"), rtCopula)

setMethod("pCopula", signature("matrix", "tCopula"), ptCopula)
setMethod("pCopula", signature("numeric", "tCopula"),ptCopula)
setMethod("dCopula", signature("matrix", "tCopula"), dtCopula)
setMethod("dCopula", signature("numeric", "tCopula"),dtCopula)

setMethod("show", signature("tCopula"), showTCopula)

setMethod("tau", "tCopula",
          function(copula) 2 * asin(copula@getRho(copula)) / pi)
setMethod("rho", "tCopula",
          function(copula) asin(copula@getRho(copula) / 2) * 6 / pi)

setMethod("lambda", signature("tCopula"), lambdaTCopula)
