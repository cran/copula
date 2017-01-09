## Copyright (C) 2016 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

##################################################################################
### Rotated copulas created from an existing copula and a mask of logicals
##################################################################################

setClass("rotCopula", contains = "xcopula",
	 slots = c(flip = "logical"), # plus "copula"
	 validity =
	     function(object) {
		 d <- dim(object@copula)
		 if((lfl <- length(object@flip)) != 1L && lfl != d)
		     "'length(flip)' must be one or match the copula dimension"
		 else if(anyNA(object@flip))
		     "'flip' contains NA/NaN"
		 else TRUE
	     })


setClass("rotExplicitCopula", contains = "rotCopula",
         slots = c(exprdist = "expression")
         ## validity = function(object) {
         ##     if (.hasSlot(object@copula, "exprdist") &&
         ##         is.language(object@copula@exprdist)) TRUE
         ##     else "The copula is not an explicit copula."
         ## }
         )



##' @title Rotated Copulas Created from an Existing Copula and a Mask of Logicals
##' @param copula The 'base' copula
##' @param flip logical (vector); if element i is TRUE,
##'        the copula is "rotated" wrt the axis x_i = 0.5;
##'        the default value is all TRUE which gives the survival copula
##' @return a new "rotCopula" object; see above
##' @author Ivan Kojadinovic (and Martin Maechler)
rotCopula <- function(copula, flip = TRUE) {
    if (isExplicit(copula))
        rotExplicitCopula(copula, flip)
    else {
        if (inherits(copula, "rotCopula")) {
            copula@flip <- copula@flip != flip
            copula
        } else
            new("rotCopula", copula = copula, flip = flip)
    }
}

##' @title Rotated Explicit Copulas Created from an Existing Copula and a Mask of Logicals
##' @param copula The 'base' copula
##' @param flip logical (vector); if element i is TRUE,
##'        the copula is "rotated" wrt the axis x_i = 0.5;
##'        the default value is all TRUE which gives the survival copula
##' @return a new "rotExplicitCopula" object; see above
##' @author Jun Yan

rotExplicitCopula <- function(copula, flip = TRUE) {
    stopifnot(isExplicit(copula))
    d <- dim(copula)
    ## TODO: if (inherits(copula, "rotExplicitCopula")) then just "flip the flips"
    if (length(flip) == 1) flip <- rep(flip, d)
    else if(length(flip) != d) stop("length(flip) must be 1 or d")

    ## preparation for cdf
    cdf <- copula@exprdist$cdf
    cdf <- do.call(substitute, list(cdf, list(alpha = quote(param))))
    unames <- paste0("u", 1L:d)
    lo <- ifelse(flip, unames, 0L)
    up <- ifelse(flip, 1L, unames)
    ## follow prob to construct the cdf expression
    D <- 2^d
    m <- 0:(D - 1)
    ## digitsBase() from package 'sfsmisc' {slightly simplified} :
    ## Purpose: Use binary representation of 0:N
    ## Author: Martin Maechler, Date:  Wed Dec  4 14:10:27 1991
    II <- matrix(0, nrow = D, ncol = d)
    for (i in d:1L) {
        II[,i] <- m %% 2L + 1L
        if (i > 1) m <- m %/% 2L
    }
    ## Sign: the ("u","u",...,"u") case has +1; = c(2,2,...,2)
    Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2]
    U <- array(cbind(lo, up)[cbind(c(col(II)), c(II))], dim = dim(II))
    ## sum(Sign * pCopula(U, x))
    rotCdf <- 0
    for (i in 1:D) {
        rep.l <- parse(text = paste0(
                           "list(",
                           paste0(unames, " = quote(", U[i, ], ")", collapse = ", "),
                           ")"))
        term <- do.call(substitute, list(cdf, eval(rep.l)))
        rotCdf <- substitute(a + sgn * b, list(a = rotCdf, b = term, sgn = Sign[i]))
    }
    oldu <- paste0("u", 1:d)[flip]
    newu <- paste0("1 - ", oldu)
    flip.l <- parse(text = paste0(
                        "list(",
                        paste0(oldu, " = quote(", newu, ")", collapse = ", "),
                        ")"))
    rotCdf <- do.call(substitute, list(rotCdf, eval(flip.l)))

    cdf <- as.expression(rotCdf)
    cdf.algr <- deriv(cdf, "nothing")

    ## preparation for pdf
    ## if (inherits(copula, "rotExplicitCopula")) {
    ##     u <- paste0("u", 1L:d)
    ##     omu <- paste0("1 - u", 1L:d) ## one minus u
    ##     oldu <- ifelse( copula@flip, omu, u)[flip]
    ##     newu <- ifelse(!copula@flip, u, omu)[flip]
    ## } else {
    ## }
    pdf <- copula@exprdist$pdf
    pdf <- do.call(substitute, list(pdf, list(alpha = quote(param))))
    pdf <- as.expression(do.call(substitute, list(pdf, eval(flip.l))))
    pdf.algr <- deriv(pdf, "nothing")

    exprdist <- c(cdf = cdf, pdf = pdf)
    attr(exprdist, "cdfalgr") <- cdf.algr
    attr(exprdist, "pdfalgr") <- pdf.algr
    new("rotExplicitCopula", copula = copula, flip = flip,
        exprdist = exprdist)
}


##################################################################################
### Basic methods
##################################################################################

## dimension setMethod("dim", signature("rotCopula"), ...) : via "xcopula"


## parameter names
setMethod("paramNames", signature("rotCopula"), function(x) paramNames(x@copula))

## get parameters
setMethod("getTheta", signature("rotCopula"),
	  function(copula, freeOnly = TRUE, attr = FALSE, named = attr)
	      getTheta(copula@copula, freeOnly=freeOnly, attr=attr, named=named))

## set free parameters
setMethod("freeParam<-", signature("rotCopula", "numeric"),
          function(copula, value) {
    freeParam(copula@copula) <- value
    copula
})

## set or modify "fixedness" of parameters
setMethod("fixedParam<-", signature("rotCopula", "logical"),
          function(copula, value) {
    fixedParam(copula@copula) <- value
    copula
})

## describe copula
setMethod(describeCop, c("rotCopula", "character"), function(x, kind, prefix="", ...)
    paste0(prefix, "Rotated copula constructed from\n", prefix, describeCop(x@copula, kind, prefix, ...)))

##################################################################################
### Methods for rotated copulas
##################################################################################

## Internal. swicth u[,i] to 1 - u[,i] according to flip
apply.flip <- function(u, flip) {
    if(identical(flip, TRUE))
        1 - u
    else if(identical(flip, FALSE))
        u
    else {
        u[,flip] <- 1 - u[,flip]
        u
    }
}

## pCopula
pRotCopula <- function(u, copula, ...) {
    apply(apply.flip(u, copula@flip), 1L, # TODO: vectorize prob ?
          function(x) prob(copula@copula,
                           l = pmin(x, copula@flip),
                           u = pmax(x, copula@flip)))
}

## dCopula
dRotCopula <- function(u, copula, log = FALSE, ...) {
    dCopula(apply.flip(u, copula@flip), copula@copula, log = log, ...)
}

setMethod("pCopula", signature("matrix",  "rotCopula"), pRotCopula)
setMethod("dCopula", signature("matrix",  "rotCopula"), dRotCopula)
setMethod("pCopula", signature("matrix",  "rotExplicitCopula"), pExplicitCopula.algr)
setMethod("dCopula", signature("matrix",  "rotExplicitCopula"), dExplicitCopula.algr)
## pCopula() and dCopula() *generic* already deal with non-matrix case!

## rCopula
setMethod("rCopula", signature("numeric", "rotCopula"),
          function(n, copula) apply.flip(rCopula(n, copula@copula), copula@flip))

##' sign() of a bivariate rotated copula
##' @return -1 : if the two 'flip' differ;  +1 : if they are (T,T) or (F,F)
sign.rot2C <- function(flip)
    if(length(flip) == 1L || flip[1L] == flip[2L]) +1 else -1

## rho
setMethod("rho", signature("rotCopula"), function(copula) {
    stopifnot(dim(copula@copula) == 2)
    sign.rot2C(copula@flip) * rho(copula@copula)
})

## iRho
setMethod("iRho", signature("rotCopula"), function(copula, rho)
    iRho(copula@copula, sign.rot2C(copula@flip) * rho))

## tau
setMethod("tau", signature("rotCopula"), function(copula) {
    sign.rot2C(copula@flip) * tau(copula@copula)
})

## iTau
setMethod("iTau", signature("rotCopula"), function(copula, tau) {
    iTau(copula@copula, sign.rot2C(copula@flip) * tau)
})

## lambda
setMethod("lambda", signature("rotCopula"), function(copula) {
    stopifnot(dim(copula@copula) == 2)
    sm <- sum(rep(copula@flip, length.out = 2))
    if (sm == 1) {
        warning("lambda() method for copula class 'rotCopula' not implemented yet")
        c(lower= NA, upper= NA)
    }
    else {
        ti <- lambda(copula@copula)
        names(ti) <- NULL
        if (sm == 0)
            c(lower = ti[1], upper = ti[2])
        else ## sm == 2
            c(lower = ti[2], upper = ti[1])
    }
})

## dCdu: Restricted to *bivariate* rotated copulas
## Needed for Khoudraji copula density evaluation
setMethod("dCdu", signature("rotCopula"),
          function(copula, u, ...) {
    if (dim(copula) > 2)
        stop("method 'dCdu' not implemented for 'rotCopula' objects in dimension higher than two")
    if (length(copula@flip) == 1L)
        apply.flip(dCdu(copula@copula, apply.flip(u, copula@flip)), copula@flip)
    else
        apply.flip(dCdu(copula@copula, apply.flip(u, copula@flip)), copula@flip[2:1])
})

## dCdtheta: Restricted to *bivariate* rotated copulas
## Needed for Khoudraji copulas
setMethod("dCdtheta", signature("rotCopula"),
          function(copula, u, ...) {
    if (dim(copula) > 2)
        stop("method 'dCdtheta' not implemented for 'rotCopula' objects in dimension higher than two")
    sign.rot2C(copula@flip) * dCdtheta(copula@copula, apply.flip(u, copula@flip))
})

## fitCopula
setMethod("fitCopula", signature("rotCopula"), function(copula, data, ...) {
    if(!is.matrix(data)) {
        warning("coercing 'data' to a matrix.")
        data <- as.matrix(data); stopifnot(is.matrix(data))
    }
    fit <- fitCopula(copula@copula, data = apply.flip(data, copula@flip), ...)
    copula@copula <- fit@copula
    fit@copula <- copula
    fit
})

## gofCopula
setMethod("gofCopula", signature("rotCopula"), function(copula, x, ...) {
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    x[,copula@flip] <- -x[,copula@flip]
    gofCopula(copula@copula, x = x, ...)
})
