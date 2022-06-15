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


### Auxiliary transformations for copula observations ##########################

##' @title Pseudo-observations
##' @param x matrix of random variates to be converted to pseudo-observations
##' @param na.last passed to rank()
##' @param ties.method passed to rank()
##' @param lower.tail if FALSE, pseudo-observations when apply the empirical
##'        marginal survival functions are returned.
##' @return pseudo-observations (of the same dimensions as x)
##' @author Marius Hofert & Martin Maechler
pobs <- function(x, na.last = "keep",
		 ## formals(rank) works in pre-2015-10-15 and newer version of rank():
		 ties.method = eval(formals(rank)$ties.method),
		 lower.tail = TRUE)
{
    ties.method <- match.arg(ties.method)
    U <- if(!is.null(dim(x)))
	     apply(x, 2, rank, na.last=na.last, ties.method=ties.method) / (nrow(x)+1)
	 else
	     rank(x, na.last=na.last, ties.method=ties.method) / (length(x)+1)
    if(inherits(x, "zoo")) # incl "xts" (but no similar) -- FIXME? and use:
### if(is.object(x) && !isS4(x) && !is.data.frame(x)) ## "zoo", "xts" et al
	attributes(U) <- attributes(x)
    if(lower.tail) U else 1-U
}

