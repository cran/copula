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


### Quantities related to the empirical copula #################################

## DEPRECATED
Cn <- function(x,w) {
    .Deprecated("C.n")
    C.n(w, x)
}

##' Empirical copula of X at x
##'
##' @title Empirical CDF and hence copula of X at x
##' @param x (m, d) matrix of evaluation points
##' @param X (n, d) matrix of pseudo-data based on which the empirical copula
##'        is computed.
##' @param offset scaling factor of result which is  sum(....)/(n+offset)
##' @param method method string ("C" for C code; "R" for R code)
##' @return empirical CDF / copula of X at x
##' @author Ivan Kojadinovic, Marius Hofert and Martin (C.n -> F.n; re-organisation)
##' Note: See ../man/empcop.Rd for a nice graphical check with the Kendall function

## Internal function
.Fn <- function(x, X, smoothing=c("none", "beta", "checkerboard"),
                offset=0, method=c("C", "R"))
## {
##     if(!is.matrix(x)) x <- rbind(x, deparse.level=0L)
##     if(!is.matrix(X)) X <- rbind(X, deparse.level=0L)
##     structure(class = "mvFn", x=x, X=X, ## <- so the result can be plotted/printed ..
##               .Fn(x, X, offset=offset, method=method))
## }

## .Fn <- function(x, X, offset=0, method=c("C", "R"))
{
    stopifnot(is.numeric(d <- ncol(X)), is.matrix(x), d == ncol(x))
    n <- nrow(X)
    if(d == 1)
	vapply(x, function(x.) sum(X <= x.), NA_real_) / (n+offset)
    else { ## d > 1
	method <- match.arg(method)
        smoothing <- match.arg(smoothing)
        type <- which(smoothing == c("none", "beta", "checkerboard"))
	switch(method,
	       "C"={
		   m <- nrow(x)
		   .C(Cn_C,		# see ../src/empcop.c
		      as.double(X),
		      as.integer(n),
		      as.integer(d),
		      as.double(x),
		      as.integer(m),
		      ec=double(m),
		      as.double(offset),
                      as.integer(type))$ec
	       },
	       "R"={
                   ## == apply(x, 1, function(x.) sum(colSums(t(X)<=x.)==d)/(n+offset) )
                   ## but vapply is slightly faster (says MH)
                   tX <- t(X)
                   switch(smoothing,
                          "none"={
                             vapply(1:nrow(x), function(k) sum(colSums(tX <= x[k,]) == d),
                                    NA_real_) / (n + offset)
                          },
                          "beta"={
                             stop("no R implementation for 'smoothing': ", smoothing)
                          },
                          "checkerboard"={
                             stop("no R implementation for 'smoothing': ", smoothing)
                          },
                          stop("wrong 'smoothing': ", smoothing)
                          )
               },
               stop("wrong 'method': ", method)
               )
    }
}

F.n <- function(x, X, offset=0, method=c("C", "R"))
    .Fn(x=x, X=X, offset=offset, method=method, smoothing="none")

C.n <- function(u, X, smoothing=c("none", "beta", "checkerboard"),
                offset=0, method = c("C", "R")) {
    if(any(u < 0, 1 < u))
        stop("'u' must be in [0,1].")
    .Fn(u, pobs(X), offset=offset, method=method, smoothing=smoothing)
}

##' @title Estimated Partial Derivatives of a Copula Given the Empirical Copula
##' @param u (m, d)-matrix of evaluation points
##' @param U (n, d)-matrix of pseudo-observations
##' @param j.ind dimensions for which the partial derivatives should be estimated
##' @param b bandwidth in (0, 1/2) for the approximation
##' @param ... additional arguments passed to F.n()
##' @return (m, length(j.ind))-matrix containing the estimated partial
##'         derivatives with index j.ind of the empirical copula of U at u
##' @author Marius Hofert
dCn <- function(u, U, j.ind=1:d, b=1/sqrt(nrow(U)), ...)
{
    ## Check
    if(!is.matrix(u)) u <- rbind(u, deparse.level=0L)
    if(!is.matrix(U)) U <- rbind(U, deparse.level=0L)
    stopifnot((d <- ncol(U)) == ncol(u),
              0 <= u, u <= 1, 0 <= U, U <= 1,
              1 <= j.ind, j.ind <= d, 0 < b, b < 0.5)

    ## Functions to change the entry in the jth column of u
    ## See Remillard, Scaillet (2009) "Testing for equality between two copulas"
    adj.u.up <- function(x){
        x. <- x + b
        x.[x > 1-b] <- 1
        x.
    }
    adj.u.low <- function(x){
        x. <- x - b
        x.[x < b] <- 0
        x.
    }

    ## (v)apply to each index...
    m <- nrow(u)
    res <- vapply(j.ind, function(j){
        ## Build the two evaluation matrices for Cn (matrix similar to u with
        ## column j replaced). Note, this is slightly inefficient for those u < b
        ## but at least we can "matricize" the problem.
        ## 1) compute F.n(u.up)
        u.up <- u
        u.up[,j] <- adj.u.up(u.up[,j])
        Cn.up <- F.n(u.up, U, ...)
        ## 2) compute F.n(u.low)
        u.low <- u
        u.low[,j] <- adj.u.low(u.low[,j])
        Cn.low <- F.n(u.low, U, ...)
        ## 3) compute difference quotient
        (Cn.up - Cn.low) / (2*b)
    }, numeric(m))
    res <- pmin(pmax(res, 0), 1) # adjust (only required for small sample sizes)
    if(length(j.ind)==1) as.vector(res) else res
}

