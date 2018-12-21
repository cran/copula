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


### Frechet--Hoeffding bound copula class ######################################

## Constructor (see ellipCopula() or archmCopula())
fhCopula <- function(family = c("upper", "lower"), dim = 2L)
{
    switch(match.arg(family), # error if invalid
           "lower" = lowfhCopula(dim = dim),
           "upper" =  upfhCopula(dim = dim))
}


### Methods ####################################################################

## describe method
setMethod(describeCop, c("fhCopula", "character"), function(x, kind, prefix="", ...) {
    kind <- match.arg(kind)
    cl <- sub("Copula$", "", class(x)) # -> "lowfh" / "upfh"
    cl <- if (cl == "lowfh") "lower Frechet-Hoeffding bound" else "upper Frechet-Hoeffding bound"
    clNam <- paste0(prefix, cl)
    if(kind == "very short") # e.g. for show() which has more parts
        return(clNam)
    d <- dim(x)
    ch <- paste(paste0(clNam, ", dim. d ="), d)
    switch(kind <- match.arg(kind),
	   short = ch,
	   long = ch,
	   stop("invalid 'kind': ", kind))
})

## dCopula method (other *Copula() methods are W/M specific)
setMethod("dCopula", signature("matrix", "fhCopula"),
	  function(u, copula, log = FALSE, ...) {
              stop("Frechet-Hoeffding bounds do not have densities")
              ## Alternatively:
    	      ## stopifnot(ncol(u) == copula@dimension)
	      ## rep.int(if(log) -Inf else 0, nrow(u)) # or even Inf or NA on the diagonal
})
