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


indepCopula <- function(dim = 2L) {
    ## get expressions of cdf and pdf
  cdfExpr <- function(n) {
    uis <- paste("u", 1:n, sep="")
    expr <- paste(uis, collapse="*")
    parse(text = expr)
  }

  pdfExpr <- function(cdf, n) {
    val <- cdf
    for (i in 1:n) {
      val <- D(val, paste("u", i, sep=""))
    }
    val
  }
  cdf <- cdfExpr((dim <- as.integer(dim)))
  pdf <- pdfExpr(cdf, dim)

  new("indepCopula",
             dimension = dim,
             exprdist = c(cdf=cdf, pdf=pdf),
             parameters = double(0),
             param.names = character(0),
             param.lowbnd = double(0),
             param.upbnd = double(0),
             message = "Independence copula")
}

AfunIndep <- function(copula, w) {
  rep(1, length(w))
}

rindepCopula <- function(copula, n) {
  dim <- copula@dimension
  matrix(runif(n * dim), nrow = n)
}

pindepCopula <- function(copula, u, log.p=FALSE) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  stopifnot (ncol(u) == copula@dimension)
  if(log.p) rowSums(log(u)) else apply(u, 1, prod)
}

dindepCopula <- function(copula, u, log=FALSE, ...) {
  if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
  stopifnot (ncol(u) == copula@dimension)
  rep.int(if(log) 0 else 1, nrow(u))
}

## kendallsTauIndepCopula <- function(copula) {
##   0
## }

## calibKendallsTauIndepCopula <- function(copula, tau) {
##   cat("No need to calibrate an independent copula.\n")
## }

## spearmansRhoIndepCopula <- function(copula) {
##   0
## }

## calibSpearmansRhoIndepCopula <- function(copula, rho) {
##   cat("No need to calibrate an independent copula.\n")
## }


setMethod("rcopula", signature("indepCopula"), rindepCopula)
setMethod("pcopula", signature("indepCopula"), pindepCopula)
setMethod("dcopula", signature("indepCopula"), dindepCopula)

setMethod("Afun", signature("indepCopula"), AfunIndep)
# setMethod("kendallsTau", signature("indepCopula"), kendallsTauIndepCopula)
# setMethod("spearmansRho", signature("indepCopula"), spearmansRhoIndepCopula)
# setMethod("tailIndex", signature("indepCopula"), tailIndexIndepCopula)

# setMethod("calibKendallsTau", signature("indepCopula"), calibKendallsTauIndepCopula)
# setMethod("calibSpearmansRho", signature("indepCopula"), calibSpearmansRhoIndepCopula)

# setMethod("tauDer", signature("indepCopula"), tauDerIndepCopula)
# setMethod("rhoDer", signature("indepCopula"), rhoDerIndepCopula)
