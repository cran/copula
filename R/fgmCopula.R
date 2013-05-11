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


## constructor #################################################################

fgmCopula <- function(param = NA_real_, dim = 2L) {
    if (!is.numeric(dim) || (dim <- as.integer(dim)) < 2)
        stop("dim should be an integer of at least 2")
    npar <- 2^dim - dim - 1
    if(is.na(param)) param <- rep(param, length.out = npar)
    else if (!is.numeric(param) && length(param) != npar)
        stop("wrong parameters")

    ## power set of {1,...,dim} in integer notation
    subsets <-  .C(k_power_set,
                   as.integer(dim),
                   as.integer(dim),
                   subsets = integer(2^dim))$subsets
    ## power set in character vector: {}, {1}, {2}, ..., {1,2}, ..., {1,...,dim}
    subsets.char <-  .C(k_power_set_char,
                        as.integer(dim),
                        as.integer(2^dim),
                        as.integer(subsets),
                        sc = character(2^dim))$sc

    ## expression of the cdf
    cdfExpr <- function(n,sc) {

        expr1 <- "u1"
        for (i in 2:n)
            expr1 <- paste(expr1, " * u", i, sep="")

        expr2 <- "1"
        for (i in (dim + 2):2^dim)
        {
            expr3 <- paste("alpha",i,sep="")
            sub <- substr(sc[i],2,nchar(sc[i])-1)
            for (j in eval(strsplit(sub,",")[[1]]))
                expr3 <- paste(expr3, " * (1 - u", j, ")", sep="")
            expr2 <- paste(expr2,"+",expr3)
        }

        expr <- paste(expr1," * (", expr2, ")")
        parse(text = expr)
    }

    ## expression of the pdf
    pdfExpr <- function(n,sc) {
        expr2 <- "1"
        for (i in (dim + 2):2^dim)
        {
            expr3 <- paste("alpha",i,sep="")
            sub <- substr(sc[i],2,nchar(sc[i])-1)
            for (j in eval(strsplit(sub,",")[[1]]))
                expr3 <- paste(expr3, " * (1 - 2 * u", j, ")", sep="")
            expr2 <- paste(expr2,"+",expr3)
        }
        parse(text = expr2)
    }

    cdf <- cdfExpr(dim,subsets.char)
    pdf <- pdfExpr(dim,subsets.char)

    ## create new object
    new("fgmCopula",
               dimension = dim,
               parameters = param,
               exprdist = c(cdf = cdf, pdf = pdf),
               param.names = paste("param",subsets.char[(dim+2):2^dim],sep=""),
               param.lowbnd = rep(-1, 2^dim - dim - 1),
               param.upbnd = rep(1, 2^dim - dim - 1),
               fullname = "Farlie-Gumbel-Morgenstern copula family")
}


### random number generation ###################################################

rfgmCopula <- function(n, copula) {
    dim <- copula@dimension
    alpha <- copula@parameters
    if (dim > 2)
        warning("random generation needs to be properly tested")
    val <- .C(rfgm,
              as.integer(dim),
	      c(rep.int(0., dim+1), alpha),
              as.integer(n),
              out = double(n * dim))$out
    matrix(val, n, dim, byrow=TRUE)
}


### cdf of the copula ##########################################################

pfgmCopula <- function(u, copula) {
    ## if(..) fails when u has NAs ; treatment now via generic pCopula()
    ## if (any(u < 0) || any(u > 1))
    ##     stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    param <- copula@parameters
    cdf <- copula@exprdist$cdf
    for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
    for (i in (dim + 2):2^dim) assign(paste("alpha", i, sep=""), param[i - dim - 1])
    eval(cdf)
}

## pdf of the copula ###########################################################

dfgmCopula <- function(u, copula, log=FALSE, ...) {
    ## if(..) fails when u has NAs ; treatment now via generic dCopula()
    ## if (any(u < 0) || any(u > 1))
    ##     stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    param <- copula@parameters
    pdf <- copula@exprdist$pdf
    for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
    for (i in (dim + 2):2^dim) assign(paste("alpha", i, sep=""), param[i - dim - 1])
    ## FIXME: improve log-case
    if(log) log(eval(pdf)) else eval(pdf)
}

## Kendall's tau

tauFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    2 * alpha / 9
}

## Spearman's rho

rhoFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    1 * alpha / 3
}

## calibration via tau

iTauFgmCopula <- function(copula, tau) {
  if (any(tau < -2/9 | tau > 2/9))
    warning("tau is out of the range [-2/9, 2/9]")
  pmax(-1, pmin(1, 9 * tau / 2))
}

## calibration via rho

iRhoFgmCopula <- function(copula, rho) {
  if (any(rho < -1/3 | rho > 1/3))
    warning("rho is out of the range [-1/3, 1/3]")
  pmax(-1, pmin(1, 3 * rho))
}


################################################################################

setMethod("rCopula", signature("numeric", "fgmCopula"), rfgmCopula)

setMethod("pCopula", signature("numeric", "fgmCopula"),pfgmCopula)
setMethod("pCopula", signature("matrix", "fgmCopula"), pfgmCopula)
setMethod("dCopula", signature("numeric", "fgmCopula"),dfgmCopula)
setMethod("dCopula", signature("matrix", "fgmCopula"), dfgmCopula)

setMethod("tau", signature("fgmCopula"), tauFgmCopula)
setMethod("rho", signature("fgmCopula"), rhoFgmCopula)
setMethod("iTau", signature("fgmCopula"), iTauFgmCopula)
setMethod("iRho", signature("fgmCopula"), iRhoFgmCopula)
