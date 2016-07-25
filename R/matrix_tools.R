## Copyright (C) 2012 Marius Hofert and Martin Maechler
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


### Tools around matrices

##' @title Construct a symmetric matrix with 1s on the diagonal from the given
##'        parameter vector
##' @param param parameter vector
##' @param d number of columns (or rows) of the output matrix
##' @return a symmetric matrix with 1s on the diagonal and the values of param
##'         filled column-wise below the diagonal (= row-wise above the diagonal)
##' @author Marius Hofert
p2P <- function(param, d = floor(1 + sqrt(2*length(param))))
{
    P <- diag(0, nrow=d)# (0 is faster to add)
    P[lower.tri(P)] <- param
    P <- P+t(P)
    diag(P) <- rep.int(1, d)
    P
}

##' @title Extract the vector of column-wise below-diagonal entries from a matrix
##' @param P matrix (typically a symmetric matrix as used for elliptical copulas)
##' @return the vector of column-wise below-diagonal entries of P (they are equal
##'         to the row-wise above-diagonal entries in case of a symmetric matrix)
##' @author Marius Hofert
##' Note: This is used "by foot" at several points in the package.
P2p <- function(P) P[lower.tri(P)]

##' @title Construct matrix Sigma from a given elliptical copula
##' @param copula copula
##' @return (d, d) matrix Sigma containing the parameter vector rho
##' @author Marius Hofert
getSigma <- function(copula)
{
    stopifnot(is(copula, "ellipCopula"))
    d <- copula@dimension
    rho <- copula@getRho(copula)
    switch(copula@dispstr,
	   "ex" = {
	       Sigma <- matrix(rho[1], nrow=d, ncol=d)
	       diag(Sigma) <- rep(1, d)
	       Sigma
	   },
	   "ar1" = {
	       rho^abs(outer(1:d, 1:d, FUN="-"))
	   },
	   "un" = {
	       p2P(rho, d)
	   },
	   "toep" = {
	       rho <- c(rho, 1)
	       ind <- outer(1:d, 1:d, FUN=function(i, j) abs(i-j))
	       diag(ind) <- length(rho)
	       matrix(rho[ind], nrow=d, ncol=d)
	   },
	   stop("invalid 'dispstr'"))
}

##' @title Find the Pairs with Smallest (or Largest) Entry in the (Lower)
##'        Triangular Area of a Symmetric Matrix
##' @param x A symmetric matrix
##' @param n Number of extreme values to be returned
##' @param method A character string indicating the method to be used
##' @param use.names A logical indicating whether colnames(x) are to be
##'        used (if not NULL)
##' @return A (n, 3)-matrix with the n largest/smallest/both
##'         values in the symmetric matrix x (3rd column) and the
##'         corresponding indices (1st and 2nd column)
##' @author Marius Hofert and Wayne Oldford
extremePairs <- function(x, n = 6, method = c("largest", "smallest", "both"),
                         use.names = FALSE)
{
    ## Checks
    if(!is.matrix(x)) x <- rbind(x, deparse.level = 0L)
    d <- ncol(x)
    method <- match.arg(method)
    stopifnot(n >= 1, d >= 2, nrow(x) == d, is.logical(use.names))

    ## Build (row, col)-matrix
    ind <- as.matrix(expand.grid(1:d, 1:d)[,2:1])
    ind <- ind[ind[,1]<ind[,2],] # pick out indices as they appear in the upper triangular matrix
    colnms <- colnames(x)
    if(use.names && !is.null(colnms))
        ind <- matrix(colnms[as.matrix(ind)], ncol = 2)

    ## Merge with entries to a (row, col, value)-data frame
    ## Note that, since x is symmetric, values of the *lower* triangular
    ## matrix as a vector matches indices in 'ind'
    val <- data.frame(ind, x[lower.tri(x)], stringsAsFactors = FALSE)

    ## Sort the data.frame according to the values and adjust the row/column names
    res <- val[order(val[,3], decreasing = TRUE),]
    colnames(res) <- c("row", "col", "value")
    rownames(res) <- NULL

    ## Now grab out the 'extreme' pairs and values
    pairs <- 1:nrow(ind) # d*(d-1)/2 pairs
    switch(method,
    "largest" = {
        stopifnot(n <= nrow(ind))
        res[head(pairs, n = n),] # from large to small
    },
    "smallest" = {
        stopifnot(n <= nrow(ind))
        res[rev(tail(pairs, n = n)),] # from small to large
    },
    "both" = {
        stopifnot(n <= floor(nrow(ind)/2))
        res[c(head(pairs, n = n), tail(pairs, n = n)),] # from large to small
    },
    stop("Wrong 'method'"))
}
