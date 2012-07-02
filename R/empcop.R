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

##' Empirical copula of x at w
##'
##' @title Empirical copula of x at w
##' @param x the data
##' @param w points where to evalute the empirical copula
##' @return the evaluations
##' @author Ivan Kojadinovic
Cn <- function(x,w) {

    p <- ncol(x)
    if (p < 2) stop("The data should be at least of dimension 2")
    if (ncol(w) != p)
      stop("The matrices 'x' and 'w' should have the same number of columns")

    n <- nrow(x)
    m <- nrow(w)

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## compute empirical copula at w
    .C(RmultCn,
       as.double(u),
       as.integer(n),
       as.integer(p),
       as.double(w),
       as.integer(m),
       ec = double(m))$ec
}
