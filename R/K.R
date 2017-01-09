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


### Kendall distribution #######################################################

## deprecated (former) Kendall distribution function K -- deprecated till 2016-06
K <- function(u, copula, d, n.MC = 0, log = FALSE){
    .Defunct("pK") # since 2016-06; deprecated before
    pK(u, copula = copula, d = d, n.MC = n.MC, log.p = log) # call the new function
}

##' @title Empirical Kendall distribution function K_{n,d} as in Lemma 1 of
##'        Genest, Neslehova, Ziegel (2011)
##' @param u evaluation points u in [0,1]
##' @param x data (in IR^d) based on which K is estimated
##' @return K_{n,d}(u)
##' @author Marius Hofert
##' Note: This is the empirical distribution function of a discrete radial part
##'       and thus K_{n,d}(0) > 0. The mass at the largest value of the support
##'       of R determines K_{n,d}(0)
Kn <- function(u, x)
{
    stopifnot(0 <= u, u <= 1, (n <- nrow(x)) >= 1, (d <- ncol(x)) >= 1)
    W <- vapply(seq_len(n), function(i) sum( colSums(t(x)<x[i,])==d ) / (n+1), NA_real_)
    ecdf(W)(u)
}

##' @title Kendall distribution function
##' @param u evaluation point(s) in [0,1]
##' @param copula acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC to evaluate the generator derivatives; otherwise the exact
##'        formula is used
##' @param log.p logical indicating whether the logarithm is returned
##' @return Kendall distribution function at u
##' @author Marius Hofert (and Martin Maechler
pK <- function(u, copula, d, n.MC = 0, log.p = FALSE)
{
    stopifnot(inherits(copula, "acopula"), 0 <= u, u <= 1)
    .pK(u, copula=copula, d=d, n.MC=n.MC, log.p=log.p)
}

## not exported :
.pK <- function(u, copula, d, n.MC = 0, log.p = FALSE)
{
    ## limiting cases
    n <- length(u)
    res <- u # hence keeping attributes, even ok for 'mpfr'
    res[is0 <- u == 0] <- if(log.p) -Inf else 0
    res[is1 <- u == 1] <- if(log.p) 0 else 1
    not01 <- seq_len(n)[!(is0 | is1)]
    n <- length(uN01 <- u[not01])
    if(n == 0) return(res)
    ## else: computations for the non-trivial  0 < uN01 < 1
    th <- copula@theta
    res[not01] <-
        if(n.MC > 0) { # Monte Carlo
            stopifnot(is.finite(n.MC))
            V <- copula@V0(n.MC, th) # vector of length n.MC
            psiI <- copula@iPsi(uN01, th) # vector of length n
                ## Former code: mean(ppois(d-1, V*psInv))
            lr <- -log(n.MC) + vapply(psiI, function(pI) lsum(cbind(ppois(d-1, V*pI, log.p=TRUE),
                                                                    deparse.level=0L)), numeric(1))
            if(log.p) lr else exp(lr)
        }
        else if(d == 1) { # d == 1
            if(log.p) log(uN01) else uN01 # K(u) = u
        } else if(d == 2) { # d == 2
            ## K(u) = u - psi^{-1}(u) / (psi^{-1})'(u)
            r <- uN01 + exp( copula@iPsi(uN01, theta=th, log=TRUE) -
                             copula@absdiPsi(uN01, th, log=TRUE) )
            if(log.p) log(r) else r
        } else { # d >= 3
            j <- seq_len(d-1)
            lpsiI. <- copula@iPsi(uN01, theta=th, log=TRUE)
            psiI. <- exp(lpsiI.)
            ## (d-1) x n  matrix
            ##  containing log( (-1)^j * psi^{(j)}(psi^{-1}(u)) ) in the j-th row
            lx <- vapply(j, function(j.) copula@absdPsi(psiI., theta=th, degree = j., log=TRUE),
                         numeric(n))
            lfac.j <- cumsum(log(j)) ## == lfactorial(j)
            lx <- (if(n==1) lx else t(lx)) + tcrossprod(j, lpsiI.) - lfac.j # (d-1) x n matrix
            ## lsum( < d x n matrix containing the logarithms of the summands of K> ) :
            ls <- lsum(rbind(log(uN01), lx)) # log(K(u))
            if(log.p) ls else pmin(exp(ls), 1) # ensure we are in [0,1] {numerical inaccuracy}

            ## NB: AMH, Clayton, Frank are numerically not quite monotone near one;
            ## --  this does not change that {but maybe slightly *more* accurate}:
            ## absdPsi. <- unlist(lapply(j, copula@absdPsi, u = psInv, theta = th,
            ##						 log = FALSE))
            ##		       sum(absdPsi.*psInv^j/factorial(j))
        } # else (d >= 3)
    res
} ## .pK()

##' @title Quantile function of the Kendall distribution function
##' @param p vector of probabilities (!)
##' @param copula acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied to evaluate K with
##'        sample size equal to n.MC; otherwise the exact formula is used
##' @param method method used for inverting K; currently:
##'        "default" : chooses a useful default
##'        "simple"  : straightforward root finding
##'        "sort"    : root finding after sorting the u's
##'        "discrete": evaluating K on u.grid, then finding approximating
##'                    quantiles based on these values
##'        "monoH.FC": evaluating K on u.grid, then approximating K via monotone
##'                    splines; finding quantiles via uniroot
##' @param u.grid default grid on which K is computed for method "discrete" and
##'        "monoH.FC"
##' @param ... additional arguments passed to uniroot() (for methods "sort",
##'        "simple", and "monoH.FC") or findInterval() (for method "discrete")
##' @return Quantile function of the Kendall distribution function at p
##' Note: - K(u) >= u  <==>  p >= K^{-1}(p)   (since  K(.) is montone)
##'       - K for smaller dimensions would also give upper bounds for K^{-1}(p),
##'         but even for d=2, there is no explicit formula for K^{-1}(p) known.
qK <- function(p, copula, d, n.MC = 0, log.p = FALSE,
               method = c("default", "simple", "sort", "discrete", "monoH.FC"),
               u.grid, ...)
{
    stopifnot(is(copula, "acopula"))
    if(log.p) stopifnot(p <= 0) else stopifnot(0 <= p, p <= 1)
    if(d==1) ## special case : K = identity
	return(p)
    ## limiting cases
    n <- length(p)
    res <- numeric(n) # all 0
    p0 <- if(log.p) -Inf else 0
    p1 <- if(log.p)   0  else 1
    res[is1 <- p == p1] <- 1
    is0 <- p == p0 # res[is0 <- p == 0] <- 0
    if(!any(not01 <- !(is0 | is1)))
	return(res)
    ## usual case:
    pN01 <- p[not01] # p's for which we have to determine quantiles
    lnot01 <- sum(not01) # may be < n
    ## computing the quantile function
    method <- match.arg(method)
    res[not01] <-
        switch(method,
               "default" =
           {
               ## Note: This is the same code as method="monoH.FC" (but with a
               ##       chosen grid)
               u.grid <- 0:128/128 # default grid
               K.u.grid <- .pK(u.grid, copula=copula, d=d, n.MC=n.MC, log.p=log.p)
               ## function for root finding
	       splFn <- splinefun(u.grid, K.u.grid, method = "monoH.FC")
	       fspl <- function(x, p) splFn(x) - p
               ## root finding
               vapply(pN01, function(p)
                      uniroot(fspl, p=p, interval=c(p0,p), ...)$root, NA_real_)

           },
               "simple" =               # straightforward root finding
           {
               ## function for root finding
               f <- function(t, p) .pK(t, copula=copula, d=d, n.MC=n.MC, log.p=log.p) - p
               ## root finding
               vapply(pN01, function(p)
                      uniroot(f, p=p, interval=c(p0,p), ...)$root, NA_real_)
           },
               "sort" =                 # root finding with first sorting p's
           {
               ## function for root finding
               f <- function(t, p) .pK(t, copula=copula, d=d, n.MC=n.MC, log.p=log.p) - p
               ## sort p's
               i.rev <- lnot01:1L # reverse, as they are typically sorted *increasingly*
               ord <- order(pN01[i.rev], decreasing=TRUE)
               pN01o <- pN01[i.rev][ord] # pN01 ordered in decreasing order
               ## deal with the first one (largest p) separately
               q <- numeric(lnot01)
               q[1] <- uniroot(f, p=pN01o[1], interval=c(p0, pN01o[1]), ...)$root # last quantile -- used in the following
               if(lnot01 > 1) {
		   ## FIXME(MM): rather use uniroot()s 'extendInt = "upX" etc ?!
                   eps <- 1e-4 # {FIXME for log.p!} ugly but due to non-monotonicity of K
                   ## otherwise: "Error... f() values at end points not of opposite sign"
                   for(i in 2:lnot01){
                       lower <- p0
                       upper <- min(pN01o[i], q[i-1]+eps)
                       if(FALSE){       # checks for debugging
                           if(lower >= upper) stop("lower=", lower, ", upper=", upper)
                           f.lower <- f(lower, pN01o[i])
                           f.upper <- f(upper, pN01o[i])
                           if(sign(f.lower*f.upper) >= 0)
                               stop("pN01o[",i,"]=", pN01o[i], ", f.lower=", f.lower, ", f.upper=", f.upper)
                       }
                       q[i] <- uniroot(f, p=pN01o[i], interval=c(lower, upper), ...)$root
                   }
               }
               ## "return" re-ordered result
               quo <- numeric(lnot01)   # vector of quantiles q
               quo[i.rev[ord]] <- q     # return result in original order
               quo
           },
               "discrete" = # evaluate K at a grid and compute approximate quantiles based on this grid
           {
	       stopifnot(0 <= u.grid, u.grid <= 1)
               K.u.grid <- .pK(u.grid, copula=copula, d=d, n.MC=n.MC, log.p=log.p)
               u.grid[findInterval(pN01, vec=K.u.grid,
                                   rightmost.closed=TRUE, ...)+1] # note: this gives quantiles according to the "typical" definition
           },
               "monoH.FC" = # root finding based on an approximation of K via monotone splines (see ?splinefun)
           {
	       stopifnot(0 <= u.grid, u.grid <= 1)
               ## evaluate K at a grid
               K.u.grid <- .pK(u.grid, copula=copula, d=d, n.MC=n.MC, log.p=log.p)
               ## function for root finding
               splFn <- splinefun(u.grid, K.u.grid, method = "monoH.FC")
               fspl <- function(x, p) splFn(x) - p
               ## root finding
               vapply(pN01, function(p)
                      uniroot(fspl, p=p, interval=c(p0,p), ...)$root, NA_real_)
           },
               stop("unsupported method ", method))
    res
}

##' @title Density of the Kendall distribution
##' @param u evaluation point(s) in (0,1)
##' @param copula acopula with specified parameter
##' @param d dimension
##' @param n.MC if > 0 a Monte Carlo approach is applied with sample size equal
##'	   to n.MC to evaluate the d-th generator derivative; otherwise the exact
##'        formula is used
##' @param log.p logical indicating whether the logarithm of the density is returned
##' @return Density of the Kendall distribution at u
##' @author Marius Hofert
dK <- function(u, copula, d, n.MC = 0, log.p = FALSE)
{
    stopifnot(is(copula, "acopula"), 0 < u, u < 1)
    th <- copula@theta
    lpsiI <- copula@iPsi(u, theta=th, log=TRUE) # log(psi^{-1}(u))
    lpsiIDabs <- copula@absdiPsi(u, theta=th, log=TRUE) # (-psi^{-1})'(u)
    ld <- lfactorial(d-1) # log((d-1)!)
    psiI <- copula@iPsi(u, theta=th) # psi^{-1}(u)
    labsdPsi <- copula@absdPsi(psiI, theta=th, degree=d, n.MC=n.MC, log=TRUE) # log((-1)^d psi^{(d)}(psi^{-1}(u)))
    res <- labsdPsi-ld+(d-1)*lpsiI+lpsiIDabs
    if(log.p) res else exp(res)
}

##' @title Random number generation for the Kendall distribution
##' @param n number of random variates to generate
##' @param copula acopula with specified parameter or onacopula
##' @param d dimension (only needed if ...)
##' @return Random numbers from the Kendall distribution
##' @author Marius Hofert
rK <- function(n, copula, d) {
    if(is(copula, c1 <- "acopula")) {
	stopifnot(d == round(d))
	copula <- onacopulaL(copula@name, list(copula@theta, 1L:d))
    } else if(!is(copula, c2 <- "outer_nacopula"))
	stop(gettextf("'copula' must be \"%s\" or \"%s\"", c1,c2), domain=NA)

    pCopula(rCopula(n, copula = copula), copula = copula)
}
