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


## Variance reduction methods for copulas

##' @title Latin hypercube sampling
##' @param u (n, d)-matrix of copula data
##' @param ... additional arguments passed to rank()
##' @return (n, d)-matrix containing the Latin Hypercube sample
##' @author Marius Hofert
rLatinHypercube <- function(u, ...)
{
    stopifnot(0 <= u, u <= 1)
    ## As pCopula(), we could use:
    ## u[] <- pmax(0, pmin(1, u))
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    n <- nrow(u)
    U <- matrix(runif(n * ncol(u)), nrow = n)
    (apply(u, 2, rank, ...) - 1 + U) / n
}

##' @title Antithetic variates
##' @param u (n, d)-matrix of copula data
##' @return (n, d, 2)-array containing u in .[,,1] and the corresponding
##'         antithetic sample in .[,,2]
##' @author Marius Hofert
rAntitheticVariates <- function(u)
{
    stopifnot(0 <= u, u <= 1)
    ## As pCopula(), we could use:
    ## u[] <- pmax(0, pmin(1, u))
    if(!is.matrix(u)) u <- rbind(u, deparse.level = 0L)
    array(c(u, 1-u), dim = c(nrow(u), ncol(u), 2))
}
