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

##################################################################################
### Virtual class of asymmetric copulas
##################################################################################

setClass("asymCopula", contains = c("copula", "VIRTUAL"))

##################################################################################
### Virtual class of asymmetric copulas constructed from two d-dimensional copulas
##################################################################################

setClass("asym2Copula", contains = c("asymCopula", "VIRTUAL"),
	 slots = c(
             copula1 = "copula",
             copula2 = "copula"
         ),
	 validity = function(object) {
	     if(object@copula1@dimension != object@copula2@dimension)
		 "The argument copulas are not of the same dimension"
	     else TRUE
         })

##################################################################################
### Potentially asymmetric d-dimensional copulas of the form C(u^{1-a}) * D(u^a)
### This construction is known as Khoudraji's device
##################################################################################

setClass("khoudrajiCopula", contains = c("asym2Copula"))
## slots *and* validity currently from 'asym2Copula' !

##################################################################################
### Bivariate Khoudraji copulas
### Can be constructed from any bivariate copulas
##################################################################################

setClass("khoudrajiBivCopula", contains = "khoudrajiCopula",
	 validity = function(object)
	     if(object@copula1@dimension != 2)
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
             exprdist = "expression",
             derExprs1 = "expression",
             derExprs2 = "expression"
         )
         ## validity = function(object) {
         ##     ## TODO: check exprdist, derExprs[12]
         ## },
         )

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
khoudrajiCopula <- function(copula1 = indepCopula(), copula2 = indepCopula(),
                            shapes = rep(NA_real_, dim(copula1))) {

    ## checks
    stopifnot((d <- dim(copula1)) == length(shapes))

    ## deal with possibly fixed attributes
    parameters <- c(copula1@parameters, copula2@parameters, shapes)
    attr(parameters, "fixed") <- c(isFixedP(copula1@parameters),
                                   isFixedP(copula2@parameters),
                                   isFixedP(shapes))

    ## if d==2, create a khoudrajiBivCopula object
    ## if d > 2
    ##          if copula1 and copula2 are explicit copulas
    ##               create a khoudrajiExplicitCopula object
    ##               (for which pdrCopula will work)
    ##          else create a khoudrajiCopula object
    ##               (for which only prCopula will work)


    ## check if copula1 and copula2 have 'exprdist' slots
    ## areNotBothExplicit <- if(is.na(match("exprdist", slotNames(copula1))) ||
    ##                       is.na(match("exprdist", slotNames(copula2))) ||
    ##                       !is.language(F1 <- copula1@exprdist$cdf) ||
    ##                       !is.language(F2 <- copula2@exprdist$cdf)) TRUE else FALSE

    ## for the moment, check if copula1 and copula2 are archmCopula, which,
    ## for the moment, implies that copula1 and copula2 have 'exprdist' slots
    areBothExplicit <- is(copula1, "archmCopula") && is(copula2, "archmCopula")

    ## d == 2 or non-explicit Khourdraji copulas
    if (d == 2 || !areBothExplicit)
        new(if (d == 2) "khoudrajiBivCopula" else "khoudrajiCopula",
            dimension = d,
            parameters = parameters,
            param.names = c(if (length(copula1@parameters) > 0)
                                paste0("c1.", copula1@param.names) else character(0),
                            if (length(copula2@parameters) > 0)
                                paste0("c2.", copula2@param.names) else character(0),
                            paste0("shape", 1:d)),
            param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, rep(0, d)),
            param.upbnd  = c(copula1@param.upbnd,  copula2@param.upbnd,  rep(1, d)),
            copula1 = copula1,
            copula2 = copula2,
            fullname = paste("Khoudraji copula constructed from: [",
                             copula1@fullname, "] and: [", copula2@fullname, "]"))
    else {

        ## Explicit Khourdraji copulas

        F1 <- copula1@exprdist$cdf
        F2 <- copula2@exprdist$cdf

        ## FIXME: not characters and parse(text=), rather expressions, substitute() ...

        ## The following block handles the cdf
        ##'' @title Replace ui with (ui^shp) for each ui in a cdf expression
        ##'' @param cdf cdf expression of a copula
        ##'' @param om logical, TRUE = use 1 - shp in place of shp
        ##'' @return a string of cdf expression with (ui^shp) or (ui^(1 - shp)) in place of ui
        getcdfchar <- function(cdf, om=FALSE) {
            ## FIXME: this only works up to dim 9; e.g., u10 could be replaced with u1^shp1
            ## -----  TODO: be smarter in gsub() --- rather do *NOT* use characters at all !!
            if (d >= 10) stop("The maximum implemented dim is 9.")
            cdf <- deparse(cdf)
            for (i in 1:d) {
                ui <- paste0("u", i)
                shpi <- paste0("shp", i)
                if (om) shpi <- paste("(1 - ", shpi, ")")
                replacement <- paste("(", ui, "^", shpi, ")")
                cdf <- gsub(ui, replacement, cdf)
            }
            cdf
        }

        ## FIXME: work with expressions F1 / F2, not chars...
        cdf1 <- gsub("alpha", "c1alpha", getcdfchar(F1, om=TRUE))
        cdf2 <- gsub("alpha", "c2alpha", getcdfchar(F2, om=FALSE))
        cdf <- parse(text = c("(", cdf1, ") * (", cdf2, ")"))
        ## cdf <- substitute((F1) * (F2), list(F1 = cdf1, F2 = cdf2))
        cdf.algr <- deriv(cdf, "nothing")
        
        ## The following block handles pdf
        ##'' @title Get pdf expression by differentiating the cdf with D iteratively
        ##'' @param cdf Expression of cdf
        ##'' @param n dimension
        ##'' @return Expression of pdf
        pdfExpr <- function(cdf, n) {
            ## This function returns the pdf expression by differentiating the cdf
            for (i in 1:n)
                cdf <- D(cdf, paste0("u", i))
            cdf
        }

        pdf <-
            if (d <= 6)
                pdfExpr(cdf, d)
            else {
                warning("The pdf is only available for dim 6 or lower.")
                NULL
            }
        pdf.algr <- deriv(pdf, "nothing")
        exprdist <- c(cdf, pdf)
        attr(exprdist, "cdfalgr") <- cdf.algr
        attr(exprdist, "pdfalgr") <- pdf.algr
        
        ## The following block handles the partial derivatives of a component copula cdf
        ## needed in the density
        ## derExprs: get the derivatives of cdf of order 1 to n
        ## That is: dC/du1, d2C/(du1 du2), d3C / (du1 du2 du3), ...
        ## WARNING: This is assuming exchangeable copula so that
        ## the ordering of arguments does not matter
        derExprs <- function(cdf, n) {
            val <- rep(as.expression(cdf), n + 1) ## the first one is cdf itself
            for (i in 1:n) {
                val[i + 1] <- as.expression(D(val[i], paste0("u", i)))
                ## FIXME: why is this as.expression necessary? Martin may know.
            }
            val
        }
        derExprs1 <- derExprs(F1, d)
        derExprs2 <- derExprs(F2, d)

        new("khoudrajiExplicitCopula",
            dimension = d,
            parameters = parameters,
            param.names = c(if (length(copula1@parameters) > 0)
                                paste0("c1.", copula1@param.names) else character(0),
                            if (length(copula2@parameters) > 0)
                                paste0("c2.", copula2@param.names) else character(0),
                            paste0("shape", 1:d)),
            param.lowbnd = c(copula1@param.lowbnd, copula2@param.lowbnd, rep(0, d)),
            param.upbnd  = c(copula1@param.upbnd,  copula2@param.upbnd,  rep(1, d)),
            copula1 = copula1,
            copula2 = copula2,
            exprdist = exprdist,
            derExprs1 = derExprs1, derExprs2 = derExprs2,
            fullname = paste("Khoudraji copula constructed from: [",
                             copula1@fullname, "] and: [", copula2@fullname, "]"))
    }
}

## In CRAN's copula up to 0.999-14 i.e  mid-2016: --> deprecated now
asymCopula <-
asymExplicitCopula <- function(shapes, copula1, copula2) {
    .Deprecated("khoudrajiCopula(c.1, c.2, shapes) -- *reordered* arguments")
    khoudrajiCopula(copula1, copula2, shapes)
}

##################################################################################
### Utility function for Khoudraji copulas
##################################################################################

## Returns shapes, copula1 and copula2 from any khoudrajiCopula object
getKhoudrajiCopulaComps <- function(object) {
    copula1 <- object@copula1
    copula2 <- object@copula2
    p1 <- length(copula1@parameters)
    p2 <- length(copula2@parameters)
    d <- object@dimension
    fixed <- attr(object@parameters, "fixed")
    shapes <- object@parameters[(p1 + p2) + 1:d]
    attr(shapes, "fixed") <- fixed[(p1 + p2) + 1:d]
    if (p1 > 0) {
        copula1@parameters <- object@parameters[1:p1]
        attr(copula1@parameters, "fixed") <- fixed[1:p1]
    }
    if (p2 > 0) {
        copula2@parameters <- object@parameters[p1 + 1:p2]
        attr(copula2@parameters, "fixed") <- fixed[p1 + 1:p2]
    }
    list(copula1 = copula1, copula2 = copula2, shapes = shapes)
}


##################################################################################
### Methods for all Khoudraji copulas
##################################################################################

## pCopula: for all Khoudraji copulas
pKhoudrajiCopula <- function(u, copula) {
    d <- copula@dimension
    tu <- if(is.matrix(u)) t(u) else matrix(u, nrow = d)
    comps <- getKhoudrajiCopulaComps(copula)
    p1 <- pCopula(t(tu^(1 - comps$shapes)), comps$copula1)
    p2 <- pCopula(t(tu^comps$shapes), comps$copula2)
    p1 * p2
}

setMethod("pCopula", signature("numeric", "khoudrajiCopula"), pKhoudrajiCopula)
setMethod("pCopula", signature("matrix", "khoudrajiCopula"),  pKhoudrajiCopula)

## rCopula: for all Khoudraji copulas
setMethod("rCopula", signature("numeric", "khoudrajiCopula"),
          function(n, copula) {
    comps <- getKhoudrajiCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    shapes <- comps$shapes
    u <- rCopula(n, copula1)
    v <- rCopula(n, copula2)
    d <- copula@dimension
    x <- matrix(NA, n, d)
    ig <- KhoudFn$ig
    for (i in seq_len(d)) {
        x[,i] <- pmax(ig(u[,i], 1 - shapes[i]), ig(v[,i], shapes[i]))
    }
    x
})

##################################################################################
### Methods for bivariate Khoudraji copulas
##################################################################################

## dCopula: Restricted to *bivariate*  copulas
dKhoudrajiBivCopula <- function(u, copula, log = FALSE, ...) {
    comps <- getKhoudrajiCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2

    ## the density can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula1)) || !hasMethod(dCdu, class(copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    part1 <- dCopula(gu1, copula1) *
        dgdu(u[,1], 1 - a1) * dgdu(u[,2], 1 - a2) *
        pCopula(gu2, copula2)
    part2 <- dC1du[,1] * dgdu(u[,1], 1 - a1) * dgdu(u[,2], a2) * dC2du[,2]
    part3 <- dC1du[,2] * dgdu(u[,2], 1 - a2) * dgdu(u[,1], a1) * dC2du[,1]
    part4 <- pCopula(gu1, copula1) * dCopula(gu2, copula2) *
        dgdu(u[,2], a2) * dgdu(u[,1], a1)

    ## FIXME: use lsum() and similar to get much better numerical accuracy for log - case
    if(log)
        log(part1 + part2 + part3 + part4)
    else    part1 + part2 + part3 + part4
}

setMethod("dCopula", signature("numeric", "khoudrajiBivCopula"), dKhoudrajiBivCopula)
setMethod("dCopula", signature("matrix", "khoudrajiBivCopula"), dKhoudrajiBivCopula)

## A: Restricted to *bivariate* Khoudraji copulas
## A: Pickands dependence function only if copula1 and copula2 are extreme-value
setMethod("A", signature("khoudrajiBivCopula"), function(copula, w) {
    comps <- getKhoudrajiCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2

    ## the A function can be computed only if the argument copulas are extreme-value copulas
    stopifnot(is(copula1, "evCopula"), is(copula2, "evCopula"))

    a1 <- comps$shapes[1];  a2 <- comps$shapes[2]
    den1 <- (1 - a1) * (1 - w) + (1 - a2) * w
    den2 <- a1 * (1 - w) + a2 * w
    t1 <- (1 - a2) * w / den1; t1 <- ifelse(is.na(t1), 1, t1)
    t2 <- a2 * w / den2; t2 <- ifelse(is.na(t2), 1, t2)
    den1 * A(copula1, t1) + den2 * A(copula2, t2)
})

## dCdu: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdu", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    comps <- getKhoudrajiCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2

    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula1)) || !hasMethod(dCdu, class(copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    g <- KhoudFn$g ; dgdu <- KhoudFn$dgdu
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula1)
    pC2gu2 <- pCopula(gu2, copula2)
    cbind(dgdu(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
          pC1gu1 * dgdu(u[,1], a1) * dC2du[,1],
          dgdu(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
          pC1gu1 * dgdu(u[,2], a2) * dC2du[,2])
})

## dCdtheta: Restricted to *bivariate* Khoudraji copulas
setMethod("dCdtheta", signature("khoudrajiBivCopula"),
          function(copula, u, ...) {
    comps <- getKhoudrajiCopulaComps(copula)
    a1 <- comps$shapes[1]
    a2 <- comps$shapes[2]
    copula1 <- comps$copula1
    copula2 <- comps$copula2

    ## dCdu can be computed only if dCdu is implemented for argument copulas
    if (!hasMethod(dCdu, class(copula1)) || !hasMethod(dCdu, class(copula2)))
        stop("The argument copulas must both have the 'dCdu()' method implemented")

    shapes <- comps$shapes
    g <- KhoudFn$g
    gu1 <- cbind(g(u[,1], 1 - a1), g(u[,2], 1 - a2))
    gu2 <- cbind(g(u[,1], a1), g(u[,2], a2))
    dC1du <- dCdu(copula1, gu1)
    dC2du <- dCdu(copula2, gu2)
    pC1gu1 <- pCopula(gu1, copula1)
    pC2gu2 <- pCopula(gu2, copula2)
    cbind(if (nFree(copula1@parameters) > 0) dCdtheta(copula1, gu1) * pC2gu2 else NULL,
          if (nFree(copula2@parameters) > 0) pC1gu1 * dCdtheta(copula2, gu2) else NULL,
          if (isFree(shapes)[1]) -log(u[,1]) * g(u[,1], 1 - a1) * dC1du[,1] * pC2gu2 +
          pC1gu1 * log(u[,1]) * g(u[,1], a1) * dC2du[,1] else NULL,
          if (isFree(shapes)[2]) -log(u[,2]) * g(u[,2], 1 - a2) * dC1du[,2] * pC2gu2 +
          pC1gu1 * log(u[,2]) * g(u[,2], a2) * dC2du[,2] else NULL)
})

##################################################################################
### dCopula method for Explicit Khoudraji copulas
##################################################################################

## This function uses the algorithmic expressions stored in the class object
dKhoudrajiExplicitCopula.algr <- function(u, copula, log=FALSE, ...) {
    dim <- copula@dimension
    stopifnot(!is.null(d <- ncol(u)), dim == d)
    comps <- getKhoudrajiCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    for (i in 1:dim) {
        assign(paste0("u", i), u[,i])
        assign(paste0("shp", i), comps$shapes[i])
    }
    ## WARNING: assuming one parameter (or no-parameter) copula components
    c1alpha <- copula1@parameters
    c2alpha <- copula2@parameters
    dens <- c(eval(attr(copula@exprdist, "pdfalgr")))
    if(log) log(dens) else dens
}

setMethod("dCopula", signature("numeric", "khoudrajiExplicitCopula"), dKhoudrajiExplicitCopula.algr)
setMethod("dCopula", signature("matrix",  "khoudrajiExplicitCopula"), dKhoudrajiExplicitCopula.algr)


#### The following is no longer used #################################################################
#### They are much slower than the one above
## > microbenchmark(dKhoudrajiExplicitCopula.algr(u, kcd3), dCopula(u, kcd3))
## Unit: milliseconds
##                                    expr       min        lq      mean    median
##  dKhoudrajiExplicitCopula.algr(u, kcd3)  2.107368  2.268947  2.660072  2.471952
##                        dCopula(u, kcd3) 60.112831 64.307831 67.056729 66.658310
##         uq      max neval
##   2.569195 19.22082   100
##  68.780618 89.72618   100


##' @title Get the power set of a sequence 1:d
##' @param d the dimension
##' @return a logical matrix of 2^d rows and d columns
##' where each row indicates one subset (TRUE indicating a component is selected)
getPowerSet <- function(d) { ## TODO: impossible for large d -- need "nextSet()" there!
    as.matrix(unname(expand.grid(replicate(d, list(c(TRUE,FALSE))),
                                 KEEP.OUT.ATTRS = FALSE)))
}

##' @title The derivatives of CDF of a copula needed in Khoudraji's copula density
##' @param idx logical vector of d dimension, TRUE means derivative requested;
##' it is designed to be one row from the returned matrix from getPowerSet.
##' @param u matrix of observations with d columns at which the derivatives are needed
##' @param dg matrix of derivatives of g at u
##' @param copula a copula object
##' @param derExprs expressions of length d + 1, the first one being the cdf of the copula,
##' the second one dC / du1, the third one dC2 / (du1 du2), and so on.
##' @return a vector of the derivatives requested via idx
##' @author Jun Yan
##' __NOT_EXPORTED__ ==> not checking arguments
densDers <- function(idx, u, dg, copula, derExprs) {
    ## assuming exchangeable copula1 and copula2
    ## IK: assuming one-parameter copulas also
    dorder <- sum(idx)
    alpha <- copula@parameters[1] # possibly needed in 'derExprs' below
    d <- copula@dimension
    newidx <- c((1:d)[idx], (1:d)[!idx])
    u <- u[, newidx]
    for (i in 1:d) assign(paste0("u", i), u[,i])
    dgu <- if (sum(idx) == 0) 1 else apply(dg[,idx,drop=FALSE], 1, prod)
    c(eval(derExprs[dorder + 1])) * dgu
}


dKhoudrajiExplicitCopula <- function(u, copula, log = FALSE, ...) {
    d <- copula@dimension
    if (!is.matrix(u)) u <- matrix(u, ncol = d)
    comps <- getKhoudrajiCopulaComps(copula)
    copula1 <- comps$copula1
    copula2 <- comps$copula2
    a <- matrix(comps$shapes, nrow(u), ncol(u), byrow=TRUE)
    u1 <- u ^ (1 - a)
    u2 <- u ^ a
    dg1 <- (1 - a) * u^(-a)
    dg2 <- a * u^(a - 1)
    powerSet <- getPowerSet(d)
    dens <- 0
    for (i in 1:nrow(powerSet)) {
        idx1 <- c(powerSet[i,])
        idx2 <- c(!powerSet[i,])
        ## WARNING: Again, this works for exchangeable copula components only
        part1 <- densDers(idx1, u1, dg1, copula1, copula@derExprs1)
        part2 <- densDers(idx2, u2, dg2, copula2, copula@derExprs2)
        dens <- dens + part1 * part2
        ## print(part1); print(part2)
    }
    if(log) log(dens) else dens
}
  


