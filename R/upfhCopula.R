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


### Upper Frechet--Hoeffding bound copula class ################################

## Constructor
upfhCopula <- function(dim = 2L)
{
    cdfExpr <- function(d) {
        uis <- paste0("u", 1:d, collapse = ",")
        expr <- paste0("min(", uis, ")")
        parse(text = expr)
    }
    new("upfhCopula",
        dimension = as.integer(dim),
        exprdist = c(cdf = cdfExpr(dim),
                     pdf = expression())) # FIXME? empty pdf disappears !?
}


### Methods ####################################################################

## for dCopula(), see fhCopula.R
setMethod("pCopula", signature("matrix", "upfhCopula"),
	  function(u, copula, log.p = FALSE) {
              stopifnot(ncol(u) == copula@dimension)
              res <- apply(u, 1, min)
	      if(log.p) log(res) else res
	  })
setMethod("rCopula", signature("numeric", "upfhCopula"),
          function(n, copula) matrix(runif(n), nrow = n, ncol = copula@dimension))

setMethod("tau", signature("upfhCopula"), function(copula) 1)
setMethod("rho", signature("upfhCopula"), function(copula) 1)
setMethod("lambda", signature("upfhCopula"), function(copula) c(lower = 1, upper = 1))
