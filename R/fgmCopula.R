#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################

#########################################################

## constructor

fgmCopula <- function(param, dim = 2) {
    if (!is.numeric(dim) || dim < 2)
        stop("dim should be a numeric greater than 2")    
    if (!is.numeric(param) && length(param) != 2^dim - dim - 1)
        stop("wrong parameters")

    ## power set of {1,...,dim} in integer notation
    subsets <-  .C("k_power_set", 
                   as.integer(dim),
                   as.integer(dim),
                   subsets = integer(2^dim),
                   PACKAGE="copula")$subsets
    ## power set in character vector: {}, {1}, {2}, ..., {1,2}, ..., {1,...,dim}
    subsets.char <-  .C("k_power_set_char", 
                        as.integer(dim),
                        as.integer(2^dim),
                        as.integer(subsets),
                        sc = character(2^dim),
                        PACKAGE="copula")$sc

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
    val <- new("fgmCopula",
               dimension = dim,
               parameters = param,
               exprdist = c(cdf = cdf, pdf = pdf),
               param.names = paste("param",subsets.char[(dim+2):2^dim],sep=""),
               param.lowbnd = rep(-1, 2^dim - dim - 1),
               param.upbnd = rep(1, 2^dim - dim - 1),
               message = "Farlie-Gumbel-Morgenstern copula family")
    val
}

#########################################################

## random number generation

rfgmCopula <- function(copula, n) {
    dim <- copula@dimension
    alpha <- copula@parameters
    if (dim > 2)
        warning("random generation needs to be properply tested")
    val <- .C("rfgm",
              as.integer(dim),
              as.double(c(rep(0,dim+1),alpha)),
              as.integer(n),
              out = double(n * dim),
              PACKAGE="copula")$out
    matrix(val, n, dim, byrow=TRUE)
}

#########################################################

## cdf of the copula

pfgmCopula <- function(copula, u) {
    if (any(u < 0) || any(u > 1))
        stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    if (is.vector(u)) u <- matrix(u, ncol = dim)
    param <- copula@parameters
    cdf <- copula@exprdist$cdf
    for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
    for (i in (dim + 2):2^dim) assign(paste("alpha", i, sep=""), param[i - dim - 1])
    val <- eval(cdf)
    val
}

#########################################################

## pdf of the copula

dfgmCopula <- function(copula, u) {
    if (any(u < 0) || any(u > 1))
        stop("u values should lie between 0 and 1")
    dim <- copula@dimension
    if (is.vector(u)) u <- matrix(u, ncol = dim)
    param <- copula@parameters
    pdf <- copula@exprdist$pdf
    for (i in 1:dim) assign(paste("u", i, sep=""), u[,i])
    for (i in (dim + 2):2^dim) assign(paste("alpha", i, sep=""), param[i - dim - 1])
    val <- eval(pdf)
    val
}

#########################################################

kendallsTauFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    2 * alpha / 9                              
}

#########################################################

spearmansRhoFgmCopula <- function(copula) {
    alpha <- copula@parameters[1]
    1 * alpha / 3                              
}

#########################################################

calibKendallsTauFgmCopula <- function(copula, tau) {
  if (any(tau < -2/9 | tau > 2/9))
    warning("tau is out of the range [-2/9, 2/9]")
  ifelse(tau <= -2/9, -1,
         ifelse(tau >= 2/9, 1, 9 * tau / 2))
}

#########################################################

calibSpearmansRhoFgmCopula <- function(copula, rho) {
  if (any(rho < -1/3 | rho > 1/3))
    warning("rho is out of the range [-1/3, 1/3]")
  ifelse(rho <= -1/3, -1,
         ifelse(rho >= 1/3, 1, 3 * rho))
}

#########################################################

setMethod("rcopula", signature("fgmCopula"), rfgmCopula)
setMethod("pcopula", signature("fgmCopula"), pfgmCopula)
setMethod("dcopula", signature("fgmCopula"), dfgmCopula)
setMethod("kendallsTau", signature("fgmCopula"), kendallsTauFgmCopula)
setMethod("spearmansRho", signature("fgmCopula"), spearmansRhoFgmCopula)
setMethod("calibKendallsTau", signature("fgmCopula"), calibKendallsTauFgmCopula)
setMethod("calibSpearmansRho", signature("fgmCopula"), calibSpearmansRhoFgmCopula)
