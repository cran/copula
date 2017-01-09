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
### Virtual class of asymmetric copulas
##################################################################################

setClass("asymCopula", contains = c("parCopula", "VIRTUAL"))

##################################################################################
### Virtual class of asymmetric copulas constructed from two d-dimensional copulas
##################################################################################

setClass("asym2Copula", contains = c("asymCopula", "VIRTUAL"),
	 slots = c(
             copula1 = "parCopula",
             copula2 = "parCopula",
             shapes  = "numeric"
         ),
	 validity = function(object) {
             if((d <- dim(object@copula1)) != dim(object@copula2))
                 "The argument copulas are not of the same dimension"
             else if(d != length(object@shapes))
                 "The length of 'shapes' is not equal to the dimension of the copulas"
             else TRUE
         })

##################################################################################
### Potentially asymmetric d-dimensional copulas of the form C(u^{1-a}) * D(u^a)
### This construction is known as Khoudraji's device
##################################################################################

setClass("khoudrajiCopula", contains = "asym2Copula")
## slots *and* validity currently from 'asym2Copula' !

##################################################################################
### Bivariate Khoudraji copulas
### Can be constructed from any bivariate copulas
##################################################################################

setClass("khoudrajiBivCopula", contains = "khoudrajiCopula",
	 validity = function(object)
             if(dim(object@copula1) != 2)
                 "The argument copulas must be of dimension two"
             else TRUE
         )

################################################################################
### Khoudraji *explicit* copulas
### Explicit refers to the fact that the two d-dimensional copulas
### have explicit cdfs and pdfs
################################################################################

setClass("khoudrajiExplicitCopula", contains = "khoudrajiCopula",
	 slots = c(
             exprdist  = "expression"
             ## ,
             ##  derExprs1 = "expression",
             ##  derExprs2 = "expression"
         )
         ## validity = function(object) {
         ##     ## TODO: check exprdist, derExprs[12]
         ## },
         )
##################################################################################
### Basic methods
##################################################################################

## dimension
setMethod("dim", signature("khoudrajiCopula"), function(x) dim(x@copula1))

## logical indicating which parameters are free
setMethod("isFree", signature("khoudrajiCopula"), function(copula)
    c(isFree(copula@copula1), isFree(copula@copula2), isFreeP(copula@shapes)))

## number of (free / all) parameters :
setMethod("nParam", signature("khoudrajiCopula"), function(copula, freeOnly=FALSE)
    nParam(copula@copula1, freeOnly=freeOnly) +
    nParam(copula@copula2, freeOnly=freeOnly) +
    (if(freeOnly) nFree else length)(copula@shapes))

## parameter names for freeOnly parameters
setMethod("paramNames", signature("khoudrajiCopula"), function(x) {
    c(if(nParam(x@copula1, freeOnly=TRUE) > 0L) paste0("c1.", paramNames(x@copula1)),
      if(nParam(x@copula2, freeOnly=TRUE) > 0L) paste0("c2.", paramNames(x@copula2)),
      paste0("shape", 1L:dim(x))[isFreeP(x@shapes)])
})

## get parameters
setMethod("getTheta", signature("khoudrajiCopula"),
          function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
    par1 <- getTheta(copula@copula1, freeOnly = freeOnly, attr = attr, named = named)
    par2 <- getTheta(copula@copula2, freeOnly = freeOnly, attr = attr, named = named)
    fixed <- attr(copula@shapes, "fixed")
    d <- dim(copula)
    ns <- if (freeOnly) nFree(copula@shapes) else d
    sel <- if (!is.null(fixed) && freeOnly) !fixed else TRUE
    par <- if (named)
	       c(if (length(par1) > 0) setNames(par1, paste0("c1.", names(par1))) else NULL,
		 if (length(par2) > 0) setNames(par2, paste0("c2.", names(par2))) else NULL,
		 setNames(copula@shapes[sel], paste0("shape", 1:d)[sel]))
	   else c(par1, par2, copula@shapes[sel])
    if (!attr)
	par
    else
	structure(par,
		  param.lowbnd = c(attr(par1, "param.lowbnd"),
				   attr(par2, "param.lowbnd"), rep(0, ns)),
		  param.upbnd  = c(attr(par1, "param.upbnd"),
				   attr(par2, "param.upbnd"), rep(1, ns)))
})

## set free parameters
setMethod("freeParam<-", signature("khoudrajiCopula", "numeric"),
          function(copula, value) {
    n1 <- nParam(copula@copula1, freeOnly=TRUE)
    n2 <- nParam(copula@copula2, freeOnly=TRUE)
    ns <- nFree(copula@shapes)
    if (n1 + n2 + ns != length(value))
        stop("the length of 'value' is not equal to the number of free parameters")
    if (n1 > 0L) freeParam(copula@copula1) <- value[1:n1]
    if (n2 > 0L) freeParam(copula@copula2) <- value[n1 + 1:n2]
    if (ns > 0L) {
        fixed <- !isFreeP(copula@shapes)
        copula@shapes[!fixed] <- value[n1 + n2 + 1:ns]
    }
    copula
})

## set parameters
setMethod("setTheta", "khoudrajiCopula",
          function(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly=TRUE, ...) {
    stopifnot(is.numeric(value) | (ina <- is.na(value)))
    if(any(ina)) {
        if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
        ## vectorized (and partial)  value <- NA_real_
        if(!is.double(value)) storage.mode(value) <- "double"
    }
    n1 <- nParam(x@copula1, freeOnly=freeOnly)
    n2 <- nParam(x@copula2, freeOnly=freeOnly)
    ns <- (if(freeOnly) nFree else length)(x@shapes)
    if (n1 + n2 + ns != length(value))
        stop(if(freeOnly)
                 "'length(value)' is not equal to the number of free parameters"
             else
                 "'length(value)' is not equal to the number of parameters")
    if (n1 > 0L) x@copula1 <- setTheta(x@copula1, value[1:n1],
                                       na.ok = na.ok, noCheck = noCheck, freeOnly = freeOnly, ...)
    if (n2 > 0L) x@copula2 <- setTheta(x@copula2, value[n1 + 1:n2],
                                       na.ok = na.ok, noCheck = noCheck, freeOnly = freeOnly, ...)
    valshapes <- value[n1 + n2 + 1:ns]
    if(all(ina.s <- is.na(valshapes)) || noCheck ||
       all(ina.s | (0 <= valshapes & valshapes <= 1)))
        ## parameter constraints are fulfilled
        x@shapes[if(freeOnly) isFreeP(x@shapes) else seq_along(valshapes)] <- valshapes
    else
        stop(gettextf("some shapes (=%s) are not between 0 and 1",
                      format(valshapes)), domain=NA)
    x
})

## set or modify "fixedness" of parameters
setMethod("fixedParam<-", signature("khoudrajiCopula", "logical"),
function(copula, value) {
    stopifnot(length(value) %in% c(1L, nParam(copula)))
    ## JY: seems not needed?
    ## if (identical(value, FALSE) || !any(value))
    ##    copula
    if (anyNA(getTheta(copula, freeOnly = FALSE)[value])) stop("Fixed parameters cannot be NA.")
    n1 <- nParam(copula@copula1)
    n2 <- nParam(copula@copula2)
    if (identical(value, FALSE) || !any(value)) {
         if (n1 > 0L) fixedParam(copula@copula1) <- FALSE
         if (n2 > 0L) fixedParam(copula@copula2) <- FALSE
         attr(copula@shapes, "fixed") <- NULL
    } else if (identical(value, TRUE) || all(value)) {
         if (n1 > 0L) fixedParam(copula@copula1) <- TRUE
         if (n2 > 0L) fixedParam(copula@copula2) <- TRUE
         attr(copula@shapes, "fixed") <- TRUE
    } else {
        if (n1 > 0L) fixedParam(copula@copula1) <- value[1:n1]
        if (n2 > 0L) fixedParam(copula@copula2) <- value[n1 + 1:n2]
        ns <- length(copula@shapes)
        attr(copula@shapes, "fixed") <-
            if (!any(v12 <- value[n1 + n2 + 1:ns])) NULL
            else if (all(v12)) TRUE else v12
    }
    copula
})

## describe copula
setMethod(describeCop, c("khoudrajiCopula", "character"), function(x, kind, prefix="", ...) {
    switch(kind <- match.arg(kind),
           "very short" = paste0(prefix, "Khoudraji copula constructed from\n",
                                 prefix, describeCop(x@copula1, "very short",
                                                     paste0(prefix, " ")),
                                 "\n", ## "\nand\n",
                                 prefix, describeCop(x@copula2, "very short",
                                                     paste0(prefix, " "))),
           "short" = paste0(prefix, "Khoudraji copula, dim. d = ", dim(x),
                            ", constructed from\n",
                            prefix, describeCop(x@copula1, "very short",
                                                paste0(prefix, " ")),
                            "\n", ## "\nand\n",
                            prefix, describeCop(x@copula2, "very short",
                                                paste0(prefix, " "))),
           "long" = paste0(prefix, "Khoudraji copula constructed from\n",
                           prefix, describeCop(x@copula1, kind, paste0(prefix, " ")),
                           "\n", ## "\nand\n",
                           prefix, describeCop(x@copula2, kind, paste0(prefix, " "))))
})

##################################################################################
### Generators for Khoudraji copulas
##################################################################################

## C(u_1^{1-a_1}, u_2^{1-a_2}) * D(u_1^a_1, u_2^a_1) = C(g(u, 1-a)) * D(g(u, a))
KhoudFn <-
    list(
        g = function(u, a) u^a,
        ## inverse of g :
        ig = function(u, a) u^(1 / a),
        ## derivative wrt u :
        dgdu = function(u, a) a * u^(a - 1)
        )

##################################################################################
### Constructor of Khoudraji copulas
##################################################################################

##' Creates a khoudrajiBivCopula object, a khoudrajiExplicitCopula
##' or a khoudrajCopula object
##'
##' @title Creates a khoudrajiBivCopula object, a khoudrajiExplicitCopula
##' or a khoudrajCopula object
##' @param copula1 a copula
##' @param copula2 a copula
##' @param shapes a numeric of length dim(copula) with elements in [0,1]
##' @return a new khoudrajiBivCopula, khoudrajiExplicitCopula
##' or a khoudrajCopula object
##' @author Jun Yan and Ivan Kojadinovic
##'
khoudrajiCopula <- function(copula1 = indepCopula(), copula2 = indepCopula(dim=d),
                            shapes = rep(NA_real_, dim(copula1))) {

    d <- dim(copula1)

    ## if both explicit, creat a khoudrajiExplicitCopula object
    ##    (for which pdrCopula will work)
    ## else if d==2, create a khoudrajiBivCopula object
    ##         (for which pdrCopula will work)
    ##      else (d > 2) create a khoudrajiCopula object
    ##         (for which only prCopula will work)


    ## check if copula1 and copula2 have 'exprdist' slots
    areBothExplicit <- isExplicit(copula1) && isExplicit(copula2)

    ## non-explicit Khourdraji copulas
    if (!areBothExplicit)
        new(if (d == 2) "khoudrajiBivCopula" else "khoudrajiCopula",
            copula1 = copula1,
            copula2 = copula2,
            shapes = shapes)
    else ## both components are explicit
        khoudrajiExplicitCopula(copula1, copula2, shapes)
}


##' @title Check if a copula is explicit
##' @param copula A copula object
##' @return TRUE if explicit FALSE otherwise
##' @author Jun Yan
isExplicit <- function(copula) {
    .hasSlot(copula, "exprdist") && is.language(copula@exprdist)
}


##' @title Prepare the Expressions in CDF of Khoudraji Copula
##' @param copula A copula object
##' @param prefix A character string "c1." or "c2." to rename the parameters
##' @param om Logical, standing for one minus. If TRUE 1 - u else u
##' @return An expression of C( u^shape) or C( (1 - u)^shape )
prepKhoudrajiCdfExpr <- function(copula, prefix, om = FALSE) {
    d <- dim(copula)
    cdf <- copula@exprdist$cdf
    ## originally, explicit copula expressions have alpha as parameter
    cdf <- do.call(substitute, list(cdf, list(alpha = quote(param))))
    oldParNames <- names(getTheta(copula, freeOnly=FALSE, named=TRUE)) # paramNames(copula)
    npar <- length(oldParNames)
    ## replace parameters
    if (npar > 0) {
        newParNames <- paste0(prefix, oldParNames)
        rep.l <- parse(text = paste0(
                           "list(",
                           paste0(oldParNames, " = quote(", newParNames, ")",
                                  collapse = ", "),
                           ")"))
        cdf <- do.call(substitute, list(cdf, eval(rep.l)))
    }
    ## replace ui with ui^shapei or ui^(1 - ^shapei)
    rep.l <- parse(text = paste0(
                       "list(",
                       paste0("u", 1:d, " = quote( (u", 1:d, ")^(",
                              if (om) "1 - shape" else "shape", 1:d, ") )",
                              collapse = ", "),
                       ")"))
    ## return cdf
    do.call(substitute, list(cdf, eval(rep.l)))
}

khoudrajiExplicitCopula <- function(copula1 = indepCopula(),
                                    copula2 = indepCopula(dim = d),
                                    shapes = rep(NA_real_, dim(copula1))) {
    d <- dim(copula1)
    cdf <- as.expression(substitute( (part1) * (part2),
        list(part1 = prepKhoudrajiCdfExpr(copula1, "c1.", TRUE),    # C_1(u)^ s
             part2 = prepKhoudrajiCdfExpr(copula2, "c2.", FALSE)))) # C_2(u)^(1-s)
    pdf <- if (d <= 6) cdfExpr2pdfExpr(cdf, d)
           else {
               warning("The pdf is only available for dim 6 or lower.")
               NULL
           }
    exprdist <- structure(c(cdf = cdf, pdf = pdf),
			  cdfalgr = deriv(cdf, "nothing"), # <-- FIXME; far from "optimal"
			  pdfalgr = deriv(pdf, "nothing")) # <--     (ditto)
    new("khoudrajiExplicitCopula",
        copula1 = copula1,
        copula2 = copula2,
        shapes = shapes,
        exprdist = exprdist)
}

## In CRAN's copula up to 0.999-14 i.e  mid-2016: --> deprecated now
asymCopula <-
asymExplicitCopula <- function(shapes, copula1, copula2) {
    .Deprecated("khoudrajiCopula(c.1, c.2, shapes) -- *reordered* arguments")
    khoudrajiCopula(copula1, copula2, shapes)
}

##################################################################################
### Methods for all Khoudraji copulas
##################################################################################

## pCopula: for all Khoudraji copulas
pKhoudrajiCopula <- function(u, copula, ...) {
    d <- dim(copula)
    tu <- if(is.matrix(u)) t(u) else matrix(u, nrow = d)
    p1 <- pCopula(t(tu^(1 - copula@shapes)), copula@copula1)
    p2 <- pCopula(t(tu^copula@shapes), copula@copula2)
    p1 * p2
}

setMethod("pCopula", signature("matrix", "khoudrajiCopula"), pKhoudrajiCopula)

## rCopula: for all Khoudraji copulas
setMethod("rCopula", signature("numeric", "khoudrajiCopula"),
          function(n, copula) {
    u <- rCopula(n, copula@copula1)
    v <- rCopula(n, copula@copula2)
    d <- dim(copula)
    x <- matrix(NA, n, d)
    ig <- KhoudFn$ig
    for (i in seq_len(d)) {
        x[,i] <- pmax(ig(u[,i], 1 - copula@shapes[i]),
                      ig(v[,i],     copula@shapes[i]))
    }
    x
})

##################################################################################
### Methods for bivariate Khoudraji copulas
##################################################################################

## dCopula: Restricted to *bivariate*  copulas
dKhoudrajiBivCopula <- function(u, copula, log = FALSE, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]

    ## the density can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1],     a1), g(u[,2],     a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    ddu1_ <- dgdu(u[,1], 1 - a1)
    ddu1  <- dgdu(u[,1],     a1)
    ddu2_ <- dgdu(u[,2], 1 - a2)
    ddu2  <- dgdu(u[,2],     a2)
    part1 <- dCopula(gu1, copula@copula1) * ddu1_ * ddu2_ *
        pCopula(gu2, copula@copula2)
    part2 <- dC1du[,1] * ddu1_ * ddu2 * dC2du[,2]
    part3 <- dC1du[,2] * ddu2_ * ddu1 * dC2du[,1]
    part4 <- pCopula(gu1, copula@copula1) * dCopula(gu2, copula@copula2) * ddu2 * ddu1

    ## FIXME: use lsum() and similar to get much better numerical accuracy for log - case
    if(log)
        log(part1 + part2 + part3 + part4)
    else    part1 + part2 + part3 + part4
}

setMethod("dCopula", signature("matrix", "khoudrajiBivCopula"),
          dKhoudrajiBivCopula)

## A: Restricted to *bivariate* Khoudraji copulas
## A: Pickands dependence function only if copula1 and copula2 are extreme-value
# setMethod("A", signature("khoudrajiBivCopula"), function(copula, w) {
setMethod("A", signature("khoudrajiCopula"), function(copula, w) {
    ## the A function can be computed only if the argument copulas are extreme-value copulas
    if(!is(copula@copula1, "evCopula") ||
       !is(copula@copula2, "evCopula"))
        stop("For Pickands A(<khoudrajiBivCop.>) both component copulas must be extreme-value ones")

    a1 <- copula@shapes[1];  a2 <- copula@shapes[2]
    den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
    den2 <- a1 * (1 - w) + a2 * w
    t1 <- (1 - a2) * w / den1; t1[is.na(t1)] <- 1
    t2 <-      a2  * w / den2; t2[is.na(t2)] <- 1
    den1 * A(copula@copula1, t1) + den2 * A(copula@copula2, t2)
})

## dCdu: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdu", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]

    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1],     a1), g(u[,2],     a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula@copula1)
    pC2gu2 <- pCopula(gu2, copula@copula2)
    cbind(dgdu(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 + pC1gu1 * dgdu(u[,1], a1) * dC2du[,1],
          dgdu(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 + pC1gu1 * dgdu(u[,2], a2) * dC2du[,2])
})

## dCdtheta: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdtheta", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    a1 <- copula@shapes[1]
    a2 <- copula@shapes[2]
    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula@copula1)) ||
        !hasMethod(dCdu, class(copula@copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")
    free <- isFreeP(copula@shapes)
    g <- KhoudFn$g
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula@copula1, gu1)
    dC2du <- dCdu(copula@copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula@copula1)
    pC2gu2 <- pCopula(gu2, copula@copula2)
    cbind(if(nParam(copula@copula1, freeOnly=TRUE) > 0) dCdtheta(copula@copula1, gu1) * pC2gu2,# else NULL
          if(nParam(copula@copula2, freeOnly=TRUE) > 0) pC1gu1 * dCdtheta(copula@copula2, gu2),#  "    "
          if (free[1]) -log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
          pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1],
          if (free[2]) -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
          pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2])
})

##################################################################################
### pCopula and dCopula method for Explicit Khoudraji copulas
##################################################################################

## cdf is used only for *testing* with dim = 2
pExplicitCopula.algr <- function(u, copula, log=FALSE, ...)
    .ExplicitCopula.algr(u, copula=copula, log=log, algoNm = "cdfalgr", ...)

dExplicitCopula.algr <- function(u, copula, log=FALSE, ...)
    .ExplicitCopula.algr(u, copula=copula, log=log, algoNm = "pdfalgr", ...)

setMethod("dCopula", signature("matrix",  "khoudrajiExplicitCopula"), dExplicitCopula.algr)

