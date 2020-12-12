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


### Model selection for copulas based on (k-)cross-validation #####################


##' @title Model selection for copulas based on (k-)cross-validation
##' @param copula object of type 'copula' representing the copula to be evaluated
##'        as model
##'        (if necessary, parameters will be used as starting values for fitCopula())
##' @param x (n, d)-matrix containing the data
##' @param k number of blocks for cross-validation
##'        if NULL, leave-one out cross-validation is performed
##' @param verbose logical indicating whether a progress bar is shown
##' @param ties.method passed to pobs
##' @return "cross validated log-likelihood"
##' @author Ivan Kojadinovic, Martin Maechler
xvCopula <- function(copula, x, k=NULL, verbose = interactive(),
                     ties.method = eval(formals(rank)$ties.method), ...)
{
    ## checks -- not too many! -- fitCopula() does check [and is generic!]
    ties.method <- match.arg(ties.method)
    if(!is.matrix(x))
    {
        warning("coercing 'x' to a matrix.")
        stopifnot(is.matrix(x <- as.matrix(x)))
    }
    stopifnot(is.numeric(x), (d <- ncol(x)) > 1, (n <- nrow(x)) > 0, dim(copula) == d)
    k <- if (is.null(k)) n else as.integer(k)
    stopifnot(k >= 2L, 
              (p <- n %/% k) >= 1L) # ideal size of blocks; blocks of size p+1 may exist

    if(k < n) ## shuffle lines of x  if 2 <= k < n
	x <- x[sample.int(n), ]

    ## setup progress bar
    if(verbose) {
	pb <- txtProgressBar(max=k, style=if(isatty(stdout())) 3 else 1)
	on.exit(close(pb)) # on exit, close progress bar
    }

    ## cross-validation
    xv <- 0
    m <- rep(p, k) # sizes of blocks initialised at p
    r <- n - k * p # remaining number of lines
    if (r > 0) m[seq_len(r)] <- p + 1L # size of first r blocks incremented if r > 0
    b <- c(0L, cumsum(m)) # 0 + ending line of each block
    v <- matrix(NA, p + 1, d) # points where copula density will be evaluated

    ## for each block
    for (i in seq_len(k)) {
        sel <- (b[i] + 1):b[i+1] # m[i] lines of current block
        ## estimate copula from all lines except those in sel
        u <- pobs((x.not.s <- x[-sel, , drop=FALSE]), ties.method = ties.method)
        copula <- fitCopula(copula, u, #method = "mpl",
                            estimate.variance=FALSE, ...)@copula
        imi <- seq_len(m[i])
        v.i   <- v[imi, , drop=FALSE] # (for efficiency)
        x.sel <- x[sel, , drop=FALSE]
        ## points where copula density will be evaluated
        nmi <- n - m[i] # == nr. of obs. in x.not.s
	for (j in seq_len(d)) {
	    ## v.i[,j] <- nmi * ecdf(x.not.s[,j])(x.sel[,j]) :
	    vals <- unique.default(xj <- sort(x.not.s[,j]))
	    v.i[,j] <- approxfun(vals, cumsum(tabulate(match(xj, vals))),
				 method = "constant", ties = "ordered",
				 yleft = 0, yright = nmi)(x.sel[,j])
	}
        ## rescale to *inside* (0,1) avoiding values {0,1} *symmetrically*:
        v.i <- (v.i + 1/2) / (nmi + 1)
        ## cross-validation for block i
        xv <- xv + mean(dCopula(v.i, copula, log = TRUE))

        if(verbose) ## update progress bar
            setTxtProgressBar(pb, i)
    }

    xv / k * n
}
