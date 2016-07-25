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
		 d <- object@copula@dimension
		 if((lfl <- length(object@flip)) != 1L && lfl != d)
		     "'length(flip)' must be one or match the copula dimension"
		 else if(anyNA(object@flip))
		     "'flip' contains NA/NaN"
		 else TRUE
	     })

##' @title Rotated Copulas Created from an Existing Copula and a Mask of Logicals
##' @param copula The 'base' copula
##' @param flip logical (vector); if element i is TRUE,
##'        the copula is "rotated" wrt the axis x_i = 0.5;
##'        the default value is all TRUE which gives the survival copula
##' @return a new "rotCopula" object; see above
##' @author Ivan Kojadinovic
rotCopula <- function(copula, flip = TRUE) {
    new("rotCopula", copula = copula, flip = flip)
}

## Internal. swicth u[,i] to 1 - u[,i] according to flip
apply.flip <- function(u, flip) {
    if(identical(flip, TRUE))
        u
    else if(identical(flip, FALSE))
        1 - u
    else {
        u[,flip] <- 1 - u[,flip]
        u
    }
}

## pCopula
pRotCopula <- function(u, copula) {
    apply(apply.flip(u, copula@flip), 1L, # TODO: vectorize prob ?
          function(x) prob(copula@copula,
                           l = pmin(x, copula@flip),
                           u = pmax(x, copula@flip)))
}

setMethod("pCopula", signature("numeric", "rotCopula"), pRotCopula)
setMethod("pCopula", signature("matrix", "rotCopula"), pRotCopula)

## dCopula
dRotCopula <- function(u, copula, log = FALSE, ...) {
    dCopula(apply.flip(u, copula@flip), copula@copula, log = log, ...)
}

setMethod("dCopula", signature("numeric", "rotCopula"), dRotCopula)
setMethod("dCopula", signature("matrix", "rotCopula"), dRotCopula)

## rCopula
setMethod("rCopula", signature("numeric", "rotCopula"),
          function(n, copula) apply.flip(rCopula(n, copula@copula), copula@flip))

##' sign() of a bivariate rotated copula
##' @return -1 : if the two 'flip' differ;  +1 : if they are (T,T) or (F,F)
sign.rot2C <- function(flip)
    if(length(flip) == 1L || flip[1L] == flip[2L]) +1 else -1

## rho
setMethod("rho", signature("rotCopula"), function(copula) {
    stopifnot(copula@copula@dimension == 2)
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
    stopifnot(copula@copula@dimension == 2)
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
