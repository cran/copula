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

##' Length of parameter vector for an 'ellipCopula' (tCopula or normCopula)
##'
##' @title Parameter length of an 'ellipCopula' [rho-part only!]
##' @param dim dimension of copula
##' @param dispstr dispersion string
##' @param df.fixed [for tCopula():] logical indicating if 'df' is fixed or a parameter
##' @return integer
##' @author Martin Maechler
## NOTE: This is related to  validRho() (== validityMethod for "ellipCopula") in ./Classes.R
npar.ellip <- function(dim, dispstr) {
    switch(dispstr, ## also checking for correct 'dispstr'
	   "ar1" =,
           "ex" = 1 ,
	   "un" = dim * (dim - 1) / 2,
	   "toep" = dim - 1,
	   ## otherwise
	   stop("'dispstr' not supported (yet)"))
}

## Lower bound for the rho_j
lowbnd.rho.ellip <- function(dim, dispstr, pdim = npar.ellip(dim, dispstr)) {
    switch(dispstr, ## also checking for correct 'dispstr'
	   "ex" = rep( - 1 / (dim-1), pdim),
	   "ar1"= rep( -1, pdim),
	   ## These bounds are not at all tight, but not simple box-constraints anyway:
	   "un" = ,
	   "toep" = rep(-1, pdim),# -> ~/R/MM/NUMERICS/toeplitz-posdef.R for some bounds
	   ## otherwise
	   return("'dispstr' not supported (yet)"))
}

ellipCopula <- function(family = c("normal", "t"), param = NA_real_, dim = 2L,
                        dispstr = "ex", df = 4, ...)
{
  switch(match.arg(family), # error if invalid
	 "normal" = normalCopula(param, dim = dim, dispstr = dispstr),
         "t" =           tCopula(param, dim = dim, dispstr = dispstr, df = df, ...))
}

##' For "toep"litz dispersion structure, allow specifying permutation of 1:dim indicating
##' after which permutation p[], the d-dim  (X[p[1]], X[p[2]], ..., X[p[d]])
##' has a toeplitz correlation/covariance structure
##' @return a string, possibly with "perm" attribute (so it is a valid <ellipCopula>@dispstr slot!)
dispstrToep <- function(perm = NULL, check = TRUE) {
    ans <- "toep"
    if(!is.null(perm)) {
        if(check)
            stopifnot(is.numeric(perm), (d <- length(perm)) >= 2,
                      perm == trunc(perm), range(perm) == c(1L, d))
        attr(ans, "perm") <- as.integer(perm)
    }
    # return
    ans
}


printEllipDisp <- function(x) {
    if(x@dimension > 2)
        cat("dispstr: ", x@dispstr,
            if(!is.null(p <- attr(x@dispstr, "perm")))
                paste0("(perm=(",paste(p, collapse=" "),"))"),
            "\n")
}

printEllipCopula <- function(x, ...) {
  printCopula(x, ...)
  printEllipDisp(x)
  invisible(x)
}

## traditionally still have  separate print()/show() methods for
## "tCopula"      -> ./tCopula.R       and
## "normalCopula" -> ./normalCopula.R

setMethod(describeCop, c("ellipCopula", "character"), function(x, kind, prefix="", ...) {
    kind <- match.arg(kind)
    cl <- sub("Copula$", "", class(x)) # -> "t" / "normal"
    if (cl == "normal") cl <- "Normal"
    clNam <- paste0(prefix, cl, if(nchar(cl) <= 2) "-" else " ", "copula")
    if(kind == "very short") # e.g. for show() which has more parts
        return(clNam)
    ## else
    d <- dim(x)
    ch <- paste(paste0(clNam, ", dim. d ="), d)
    switch(kind <- match.arg(kind),
	   short = ch,
	   long = paste0(ch, "\n", prefix, "param.: ", dputNamed(getTheta(x, named=TRUE))),
	   stop("invalid 'kind': ", kind))
})


iTauEllipCopula <- function(copula, tau) sinpi(tau / 2)

## iRho --> only for normalCopula

dTauEllipCopula <- function(copula)  {
  2 / (pi * sqrt(1 - copula@getRho(copula)^2))
}

dTauFunEllipCopula <- function(copula)  {
  function(x) 2 / (pi * sqrt(1 - x^2))
}

dRhoEllipCopula <- function(copula) {
  6 / (pi * sqrt(4 - copula@getRho(copula)^2))
}

dRhoFunEllipCopula <- function(copula) {
  function(x) 6 / (pi * sqrt(4 - x^2))
}

setMethod("iTau", signature("ellipCopula"), iTauEllipCopula)
## iRho: only for --> ./normalCopula.R

setMethod("dTau", signature("ellipCopula"), dTauEllipCopula)
setMethod("dRho", signature("ellipCopula"), dRhoEllipCopula)

setMethod("dTauFun", signature("ellipCopula"), dTauFunEllipCopula)
setMethod("dRhoFun", signature("ellipCopula"), dRhoFunEllipCopula)
