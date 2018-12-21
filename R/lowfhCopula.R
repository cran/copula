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


### Lower Frechet--Hoeffding bound copula class ################################

## Constructor
lowfhCopula <- function(dim = 2L)
{
    if(dim != 2) stop("The lower Frechet-Hoeffding bound copula is only available in the bivariate case")
    cdfExpr <- function(d) {
        uis <- c(paste0("u", 1:d, "-1"), "1")
        expr <- paste(uis, collapse = "+")
        expr <- paste0("max(", expr, ",0)")
        parse(text = expr)
    }
    new("lowfhCopula",
        dimension = as.integer(dim),
        exprdist = c(cdf = cdfExpr(dim),
                     pdf = expression())) # FIXME? empty pdf disappears !?
}


### Methods ####################################################################

## for dCopula(), see fhCopula.R
setMethod("pCopula", signature("matrix", "lowfhCopula"),
	  function(u, copula, log.p = FALSE) {
              d <- ncol(u)
              stopifnot(d == copula@dimension)
              res <- pmax(rowSums(u) - d + 1, 0)
	      if(log.p) log(res) else res
})
setMethod("rCopula", signature("numeric", "lowfhCopula"),
          function(n, copula) {
              U <- runif(n)
              cbind(U, 1-U)
          })

setMethod("tau", signature("lowfhCopula"), function(copula) -1)
setMethod("rho", signature("lowfhCopula"), function(copula) -1)
setMethod("lambda", signature("lowfhCopula"), function(copula) c(lower = 0, upper = 0))
