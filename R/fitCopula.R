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


### Class and methods for fitCopula ############################################

setClass("fitCopula", contains = c("xcopula", "fittedMV")) #-> ./Classes.R
##-> coef(), nobs(), vcov() and logLik() methods for 'fittedMV' there

## FIXME: This has much in commong with print.fitMvdc [ ./fitMvdc.R ]
## -----  For consistency, use common "helper" functions (instead of now: manually sync'ing)
print.fitCopula <- function(x, digits = max(3, getOption("digits") - 3),
                            signif.stars = getOption("show.signif.stars"), ...,
                            showMore = FALSE)
{
    cat("Call: ", formatCall(x@call, class(x)), "\n", sep = "")
    cop <- x@copula
    d <- dim(cop) # not slot; e.g. for rotCopula
    cat(sprintf(
	"Fit based on \"%s\" and %d %d-dimensional observations.\n",
	x@method, x@nsample, d))
    ## FIXME show more via printCopula() utility; but do *not* show the parameters
    cat(if(showMore) describeCop(cop, "short") # as 'parameters' are separate
	else paste("Copula:", class(cop)), "\n")
    if(showMore) {
	coefs <- coef.fittedMV(x, SE = TRUE)
	printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
		     na.print = "NA", ...)
    } else { # !showMore, default print() etc
	coefs <- coef.fittedMV(x) # vector of "hat(theta)"
	print(coefs, digits=digits, ...)
    }
    if (!is.na(ll <- x@loglik))
	cat("The maximized loglikelihood is", format(ll, digits=digits), "\n")
    f.stats <- x@fitting.stats
    if (!is.na(conv <- f.stats[["convergence"]])) {
	if(conv)
            cat("Convergence problems: code is", conv, "see ?optim.\n")
	else cat("Optimization converged\n")
    }
    if(showMore && !is.null(cnts <- f.stats$counts) && !all(is.na(cnts))) {
	cat("Number of loglikelihood evaluations:\n"); print(cnts, ...)
    }
    ## if(showMore)
    ##     print(cop, digits=digits, ...)
    invisible(x)
}

## NB: coef.fittedMV() does apply to <fitCopula>

## Keeping this back compatible... --> must return list(method=, loglik= ..):
summary.fitCopula <- function(object, ...) {
    structure(class = "summary.fitCopula",
              list(fitC = object,
                   method = object@method,
                   loglik = object@loglik,
                   convergence = object@fitting.stats[["convergence"]],
                   coefficients = coef.fittedMV(object, SE = TRUE)))
}

## Used both "as"  'print.summary.fitCopula()' and 'print.summary.fitMvdc()'
## Now print(summary(<fitCopula>)) *does* show a bit more than print(<fitCopula>)
printSummary.fittedMV <- function(x, ...) {
    print(x$fitC, ..., showMore = TRUE)
    invisible(x)
}

setMethod("show", signature("fitCopula"),
	  function(object) print.fitCopula(object))


### Auxiliary functions ########################################################

##' @title Construct Parameter Matrix from Matrix of
##'        Kendall's Taus / Spearman's Rhos
##' @param cop Copula object (typically with NA parameters)
##' @param x The data in [0,1]^d or IR^d
##' @param method The rank correlation measure used
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param ... Additional arguments passed to cor() or corKendall()
##' @return (d,d)-matrix of parameters *or* a d(d-1)/2 parameter vector
fitCor <- function(cop, x, method = c("itau", "irho"),
                   posDef = is(cop, "ellipCopula"), matrix = TRUE, ...) {
    method <- match.arg(method)
    theta <-
        switch(method,
        "itau" = {
            tau <- corKendall(x, ...)
            iTau(cop, P2p(tau))
        },
        "irho" = {
            rho <- cor(x, method="spearman", ...)
            iRho(cop, P2p(rho))
        },
        stop("Not yet implemented:", method))
    if (posDef) { # make pos. definite
        m <- nearPD(p2P(theta, d=ncol(x)), corr=TRUE)$mat # "dpoMatrix"
        if(!matrix) P2p(m) else m
    } else
        if(matrix) p2P(theta, d=ncol(x)) else theta
}

##' @title Computing the Log-Likelihood of the Given Copula
##' @param param Parameter (vector)
##' @param u The data in [0,1]^d
##' @param copula Copula object
##' @return Log-likelihood of the given copula at param given the data x
loglikCopula <- function(param, u, copula) {
    r <- tryCatch(freeParam(copula) <- param, error = function(e) NULL)
    if(is.null(r)) # param was presumably "out-of-bound" (TODO? look at error, maybe warn?)
        return(-Inf)
    cop.param <- getTheta(copula, freeOnly = TRUE, attr = TRUE)
    lower <- attr(cop.param, "param.lowbnd")
    upper <- attr(cop.param, "param.upbnd")
    admissible <- !any(is.na(cop.param) | cop.param > upper | cop.param < lower)
    if (admissible) {
        ## FIXME-JY: param range check is only a *part* of validity check
        sum(dCopula(u, copula=copula, log=TRUE, checkPar=FALSE))
    } else -Inf
}

##' @title Computing Initial Parameter Values for fitCopula.ml()
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d
##' @param default The default initial values for _free_ parameters
##' @param ... Additional arguments passed to fitCopula.icor()
##' @return Initial parameter value for fitCopula.ml()
fitCopStart <- function(copula, u, default = getTheta(copula, freeOnly = TRUE), ...)
{
    clc <- class(copula)
    ## start <-
    if(hasMethod("iTau", clc)) {
	ccl <- getClass(clc)
	.par.df <- has.par.df(copula, ccl)
	start <- fitCopula.icor(if(.par.df) as.df.fixed(copula, classDef = ccl) else copula,
				x=u, method="itau", estimate.variance=FALSE,
                                warn.df=FALSE, ...)@estimate # fitCopula.icor(, method="itau")
	if(.par.df) # add starting value for 'df'
	    start <- c(start, getdf(copula))
	if(is.finite(loglikCopula(start, u=u, copula=copula)))
            start
        else {
            if (is(copula, "claytonCopula") && dim(copula) == 2) {
                ## The support of bivariate claytonCopula with negative parameter is not
                ## the full unit square; the far from 0, the more restricted.
                while (start < .0) {
                    start <- start + .2
                    if (is.finite(loglikCopula(start, u=u, copula=copula))) break
                }
                start
            }
            else
                default
        }
    } else default
}

##' @title Generate the Contrast Matrix L for Sample Tau or Rho
##' @details For elliptical copulas, pairwise sample tau or rho can be pooled to
##'          estimate the true tau or rho depending on the dispersion structure.
##'          The contrast matrix L is such that the pooled estimate is expressed
##'          as t(L) multiplied by the vector of pairwise sample values. This is
##'          used in the variance estimation; see var.icor below.
##' @param copula A (most likely elliptical) copula object
##' @return L
getL <- function(copula) {
    ## for ellipCopula only!
    dim <- dim(copula) # GETR
    dd <- dim * (dim - 1) / 2
    free <- isFree(copula) # GETR
    if (.hasSlot(copula, "df.fixed")) free <- free[-length(free)]

    dgidx <- outer(1:dim, 1:dim, "-")
    dgidx <- P2p(dgidx)

    if (!is(copula, "ellipCopula") || copula@dispstr == "ex") {
        cbind(rep.int(1/dd, dd), deparse.level=0L) # no free adjustment for scalar parameter
    } else
        switch(copula@dispstr,
               "un" = diag(dd)[,free,drop=FALSE],
               "toep" = {
                   mat <- model.matrix(~ factor(dgidx) - 1)
                   mat <- mat / matrix(colSums(mat), nrow = dd, ncol= dim - 1, byrow=TRUE)
                   mat[,free,drop=FALSE]
               },
               "ar1" = {
                   ## FIXME More efficient: estimate log(rho) first and then exponentiate back,
                   ##  see e.g. fitCor(*, "irho") above
                   ## JY: can sometimes rho be negative? May need fix.
                   ## mat <- model.matrix(~ factor(dgidx) - 1)
                   ## mat * matrix(1:(dim - 1), nrow=dd, ncol=dim - 1, byrow=TRUE)
                   X <- getXmat(copula) # no free adjustment for scalar parameter
                   ## L:
                   t(solve(crossprod(X), t(X)))
               },
               stop("Not implemented yet for the dispersion structure ", copula@dispstr))
}

##' @title Design Matrix for Method-of-Moments Estimators via lm()
##' @param copula Copula object
##' @return dd by p Design matrix for method-of-moments estimators via lm()
getXmat <- function(copula) { ## FIXME(?) only works for "copula" objects, but not rotCopula, mixed*, ...
    dim <- dim(copula) # GETR
    dd <- dim * (dim - 1) / 2
    xmat <-
    if (!is(copula, "ellipCopula")) # one-parameter non-elliptical copula
	if((n.th <- nParam(copula, freeOnly=TRUE)) == 1L) # GETR
	    matrix(1, nrow=dd, ncol=1)
	else stop(gettextf( ## should not happen currently, but s
		 "getXmat() not yet implemented for non-elliptical copulas with %d parameters",
		 n.th), domain=NA)
    else { ## ellipCopula :
	switch(copula@dispstr,
	       "ex" = matrix(1, nrow=dd, ncol=1),
	       "un" = diag(dd),
	       "toep" =, "ar1" = {
                   dgidx <- P2p(outer(1:dim, 1:dim, "-"))
                   if(copula@dispstr == "toep")
                       model.matrix(~ factor(dgidx) - 1)
                   else { ## __"ar1"__
                       ## estimate log(rho) first and then exponentiate back
                       ## mat <- model.matrix(~ factor(dgidx) - 1)
                       ## mat %*% diag(1:(dim - 1))
                       cbind(dgidx, deparse.level=0L)
                   }
               },
               stop("Not implemented yet for the dispersion structure ", copula@dispstr))
    }
    free <- isFree(copula) # GETR ## the last one is for df and not needed
    if (.hasSlot(copula, "df.fixed")) free <- free[-length(free)]
    xmat[, free, drop=FALSE]
}


### Variances of the estimators ################################################

##' @title Variance of the Inversion of a Rank Correlation Measure Estimator
##' @param cop The *fitted* copula object
##' @param u The data in [0,1]^d
##' @param method A character string indicating whether Spearman's rho
##'        or Kendall's tau shall be used
##' @return Variance of the inversion of a rank correlation measure estimator
##' @note See Kojadinovic & Yan (2010) "Comparison of three semiparametric ...",
##'       IME 47, 52--63
var.icor <- function(cop, u, method=c("itau", "irho"))
{
    ## Check if variance can be computed
    method <- match.arg(method)
    dim <- dim(cop) # GETR
    dd <- dim * (dim - 1) / 2
    free <- isFree(cop) # GETR
    n <- nrow(u)
    v <- matrix(0, n, dd)

    ## Compute influence functions for the respective rank correlation measure
    if(method=="itau") {
        ec <- numeric(n)
        l <- 0L
        for (j in 1:(dim-1)) {
            for (i in (j+1):dim) {
                for (k in 1:n) # can this be vectorized? FIXME(MM)
                    ec[k] <- sum(u[,i] <= u[k,i] & u[,j] <= u[k,j]) / n
                v[,(l <- l + 1L)] <- 2 * ec - u[,i] - u[,j]
            }
        }
    } else { # Spearman's rho
        ord <- apply(u, 2, order, decreasing=TRUE)
        ordb <- apply(u, 2, rank) # ties : "average"
        storage.mode(ordb) <- "integer" # as used below
        l <- 0L
        for (j in 1:(dim-1)) {
            for (i in (j+ 1L):dim)
                v[,(l <- l + 1L)] <- u[,i] * u[,j] +
                    c(0, cumsum(u[ord[,i], j]))[n + 1L - ordb[,i]] / n +
                    c(0, cumsum(u[ord[,j], i]))[n + 1L - ordb[,j]] / n
        }
    }

    ## X <- getXmat(cop)
    ## L <- t(solve(crossprod(X), t(X)))
    L <- getL(cop)
    v <- if (is(cop, "ellipCopula") && cop@dispstr == "ar1") {
        ## Estimate log(r) first, then log(theta), and then exponentiate back
        ## r is the lower.tri of sigma
        sigma <- getSigma(cop) # assuming cop is the fitted copula
        ## Influence function for log(r)
        r <- P2p(sigma)
        D <- diag(x = 1 / r / if(method=="itau") dTauFun(cop)(r) else dRhoFun(cop)(r), dd)
        v <- v %*% D
        ## Influence function for log(theta):                 IF = v %*% L  ==>
        ## Influence function for theta = cop@parameters[1] : IF = v %*% L %*% theta
        v %*% L %*% cop@parameters[1]
    } else {
        ## Check if variance can be computed
        ## FIXME: "string" version of 'dCor' should not be needed !
        dCor <- switch(method,
                       "itau" = "dTau",
                       "irho" = "dRho")
        if (!hasMethod(dCor, class(cop))) {
            warning("The variance estimate cannot be computed for a copula of class ", class(cop))
	    return(matrix(NA_real_, 0, 0)) # instead of potentiall large (n x dd)
        }

        dCor <- if(method=="itau") dTau else dRho
        ## D <- if (length(cop@parameters) == 1) 1 / dCor(cop) else diag(1 / dCor(cop))
        D <- diag(1 / dCor(cop)[free], sum(free)) # caution: diag(0.5) is *not* a 1x1 matrix
        v %*% L %*% D
    }
    ## FIXME: can be made more efficient by extracting L and D
    ## free <- isFree(cop)
    ## v <- v[free, free, drop = FALSE]
    var(v) * if(method=="itau") 16 else 144
}

##' @title Variance-Covariance Matrix (vcov) for Pseudo Likelihood Estimate
##' @param cop The *fitted* copula object
##' @param u The data in [0,1]^d
##' @return vcov matrix
var.mpl <- function(cop, u)
{
    ## Checks
    p <- nParam(cop, freeOnly=TRUE) # parameter space dimension p
    dim <- dim(cop) # copula dimension d
    ans <- matrix(NA_real_, p, p) # matrix(NA_real_, 0, 0)
    ccl <- getClass(clc <- class(cop))
    isEll <- extends(ccl, "ellipCopula")
    ## Check if variance can be computed
    msg <- gettext("The variance estimate cannot be computed for this copula.",
                   " Rather use 'estimate.variance = FALSE'")
    if(is(cop, "archmCopula")) {
	fam <- names(which(.ac.classNames == class(cop)[[1]]))
	msg <- c(msg, gettext(" Or rather  emle(u, oCop)  instead; where",
                              sprintf(" oCop <- onacopula(%s, C(NA, 1:%d))", fam, dim)))
    }
    if (!isEll && !hasMethod("dlogcdu", clc)) {
	warning(msg); return(ans)
    }

    ## If df.fixed = FALSE, Jscore() cannot be computed
    if(has.par.df(cop, ccl, isEll)) {
        cop <- as.df.fixed(cop, classDef = ccl)
        ans[-p, -p] <- var(t(Jscore(cop, u, method = "mpl")))
        ans
    } else
        var(t(Jscore(cop, u, method = "mpl")))
}


### Estimators #################################################################

##' @title Inversion of Spearman's rho or Kendall's tau Estimator
##' @param copula The copula to be fitted
##' @param x The data in [0,1]^d or IR^d
##' @param estimate.variance A logical indicating whether the variance
##'        of the estimator shall be computed
##' @param warn.df A logical indicating whether a warning is given
##'        if the copula is coerced to df.fixed=TRUE
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param method A character string indicating whether Spearman's rho
##'        or Kendall's tau shall be used
##' @param ... Additional arguments passed to fitCor()
##' @return The fitted copula object
fitCopula.icor <- function(copula, x, estimate.variance, method=c("itau", "irho"),
                           warn.df=TRUE, posDef=is(copula, "ellipCopula"), call, ...)
{
    method <- match.arg(method)
    ccl <- getClass(class(copula))
    isEll <- extends(ccl, "ellipCopula")
    if(has.par.df(copula, ccl, isEll)) { # must treat it as "df.fixed=TRUE"
        if(warn.df)
            warning("\"", method, "\" fitting ==> copula coerced to 'df.fixed=TRUE'")
        copula <- as.df.fixed(copula, classDef = ccl)
    }
    stopifnot(any(free <- isFree(copula))) # GETR
    if(missing(call)) call <- match.call()
    if (.hasSlot(copula, "df.fixed")) free <- free[-length(free)]
    icor <- fitCor(copula, x, method=method, posDef=posDef, matrix=FALSE, ...)

    ## FIXME: Using 'X' & 'lm(X,y)' is computationally very inefficient for large q
    ## Note: The following code (possibly .lm.fit() but at least matrix(NA_real_, q, q))
    ##       can lead to R being killed (under Mac OS X) for d ~= 450
    ##       Also, it's slow in higher dimensions (even getXmat())
    estimate <- if(isEll && copula@dispstr == "un") {
        icor[free]
    } else if(isEll && copula@dispstr == "ar1") {
        X <- getXmat(copula)
        exp(.lm.fit(X, y=log(icor))$coefficients) # FIXME: assuming icor > 0
    } else {
        X <- getXmat(copula)
        .lm.fit(X, y=icor)$coefficients
    }

    estimate <- as.vector(estimate) # strip attributes
    freeParam(copula) <- estimate
    var.est <- if (is.na(estimate.variance) || estimate.variance) {
        var.icor(copula, x, method=method)/nrow(x)
    } else matrix(NA_real_, 0, 0)
    new("fitCopula",
        estimate = estimate,
        var.est = var.est,
        method = paste("inversion of", if(method=="itau") "Kendall's tau" else "Spearman's rho"),
        loglik = NA_real_,
        fitting.stats = list(convergence = NA_integer_),
        nsample = nrow(x), call = call,
        copula = copula)
}

##' @title Estimator of Mashal, Zeevi (2002) for t Copulas; see also Demarta, McNeil (2005)
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d (this would not be required if we applied pobs();
##'        the latter is fine for estimating P via pairwise inversion of Kendall's tau,
##'        but if we want a more true (up to the estimation of P) estimation of nu
##'        based on simulated copula data, this would not be possible => require the
##'        right data as input already)
##' @param posDef A logical indicating whether a proper correlation matrix
##'        is computed
##' @param lower The lower bound for optimize() (default 0 means Gaussian as we go in 1/nu)
##' @param upper The upper bound for optimize() (default 32 means down to nu = 1/32)
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed (TODO: not fully implemented yet)
##' @param tol Tolerance of optimize() for estimating nu
##' @param ... Additional arguments passed to the underlying fitCor
##' @return The fitted copula object
##' @author Marius Hofert and Martin Maechler
##' @note One could extend this to fitCopula.icor.ml(, method=c("itau", "irho"))
##'       once an *explicit* formula for Spearman's rho is available for t copulas.
fitCopula.itau.mpl <- function(copula, u, posDef=TRUE, lower=NULL, upper=NULL,
                               estimate.variance, tol=.Machine$double.eps^0.25,
                               traceOpt = FALSE, call, ...)
{
    stopifnot(any(free <- isFree(copula))) # GETR
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    ## Note: We require dispstr = "un" as it is a) the method of Mashal, Zeevi (2002) and
    ##       b) otherwise would require the use of .lm.fit which can crash for large d!
    ## MM: The above seems wrong conceptually, "un" being ill-posed for large d
    if(!(is(copula, "tCopula") && copula@dispstr == "un"))
        stop("method \"itau.mpl\" is only applicable for \"tCopula\" with 'dispstr=\"un\"'")
    if(copula@df.fixed) stop("Use method=\"itau\" for 'tCopula' with 'df.fixed=TRUE'")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (dim(copula) != d) # GETR
        stop("The dimension of the data and copula do not match")
    if(missing(call)) call <- match.call()
    if(is.null(lower)) lower <- 0  # <=>  df=Inf <=> Gaussian
    if(is.null(upper)) upper <- 32 # down to df = 1/32

    ## Estimate the correlation matrix P
    ## Note that in
    ##       fitCopula.icor(copula, u, estimate.variance=FALSE, method="itau",
    ##                      warn.df=FALSE, posDef=posDef, ...)
    ## the variance estimation fails in higher dimensions (matrix too large)
    ## matrix P ("Rho") as vector:
    P <- fitCor(copula, x = u, method = "itau", posDef = posDef, matrix=FALSE, ...)

    ## Estimate the d.o.f. parameter nu via Maximum Likelihood
    logL <-
        if(traceOpt) { ## for debugging
	    function(Inu) {
		r <- loglikCopula(c(P, df=1/Inu)[free], u=u, copula=copula) # IK: [free]
		cat(sprintf("1/nu=%14.9g => nu=%9.4f; logL=%12.8g\n",
			    Inu, 1/Inu, r))
		r
	    }
	} else
	    function(Inu) loglikCopula(c(P, df=1/Inu)[free], u=u, copula=copula) # IK: [free]

    fit <- optimize(logL, interval=c(lower, upper), tol=tol, maximum = TRUE)

    ## Extract the fitted parameters
    df <- 1/fit$maximum # '1/.' because of logL() being in '1/.'-scale
    estimate <- c(P, df)[free]
    freeParam(copula) <- estimate

    loglik <- fit$objective
    ## has.conv <- TRUE # FIXME? use tryCatch() above to catch non-convergence

    ## Deal with the variance
    if (is.na(estimate.variance))
        estimate.variance <- FALSE # not yet: using 'has.conv'
    if(estimate.variance)
        ## TODO: if (d <= 20) or so ... or  if( q < n / 5 ) then
        ## use one/zero step of fitCopula.ml(*...., method="mpl", maxit=0) to get full vcov()
        stop("Cannot estimate var-cov matrix currently for method \"itau.mpl\"")
    var.est <- matrix(NA_real_, 0, 0) # length 0  <==> not-estimated / not-available

    ## Return
    new("fitCopula",
        estimate = estimate,
        var.est = var.est,
        method = "itau for dispersion matrix P and maximum likelihood for df",
        loglik = loglik,
        fitting.stats = list(convergence = 0),
        ## optimize() does not give any info! -- if we had final optim(*, hessian=TRUE) step?
        ## c(list(method=method),
        ## fit[c("convergence", "counts", "message")], control),
        nsample = nrow(u), call = call,
        copula = copula)
}

##' @title Maximum Likelihood Estimator for Copulas
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d (for method="ml", this needs to be true copula data;
##'        for method="mpl", this can be parametrically or non-parametrically estimated
##'        pseudo-observations)
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param bound.eps A small quantity denoting an eps for staying away from
##'        the theoretical parameter bounds
##' @param ... Additional arguments (currently with no effect)
##' @return The fitted copula object
fitCopula.ml <- function(copula, u, method=c("mpl", "ml"), start, lower, upper,
                         optim.method, optim.control, estimate.variance,
                         bound.eps=.Machine$double.eps^0.5, call, traceOpt = FALSE,
			 need.finite = any(optim.method == c("Brent", "L-BFGS-B", "BFGS")),
			 finiteLARGE = .Machine$double.xmax,
			 ...)
{
    stopifnot((q <- nParam(copula, freeOnly=TRUE)) > 0L, # GETR
	      is.list(optim.control) || is.null(optim.control),
	      is.character(optim.method), length(optim.method) == 1,
              is.finite(finiteLARGE), length(finiteLARGE) == 1, finiteLARGE > 0)
    chk.s(...) # 'check dots'
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (dim(copula) != d) # GETR
        stop("The dimension of the data and copula do not match")
    if(missing(call)) call <- match.call()
    if(is.null(start))
        start <- fitCopStart(copula, u=u)
    if(anyNA(start)) stop("'start' contains NA values")
    if (q != length(start))
        stop(gettextf("The lengths of 'start' (= %d) and the number of free copula parameters (=%d) differ",
                      length(start), q), domain=NA)
    method <- match.arg(method)

    ## Determine optim() inputs
    control <- c(as.list(optim.control), fnscale = -1) # fnscale < 0 => maximization
    ## unneeded, possibly wrong (why not keep a 'special = NULL' ?):
    ## control <- control[!vapply(control, is.null, NA)]
    meth.has.bounds <- optim.method %in% c("Brent","L-BFGS-B")
    if(meth.has.bounds || need.finite)
	asFinite <- function(x) { # as finite - vectorized in 'x';  NA/NaN basically unchanged
	    if(any(nifi <- !is.finite(x)))
		x[nifi] <- sign(x[nifi]) * finiteLARGE
	    x
	}

    if(meth.has.bounds) {
        cop.param <- getTheta(copula, freeOnly = TRUE, attr = TRUE)
        p.lowbnd <- attr(cop.param, "param.lowbnd")
        p.upbnd  <- attr(cop.param, "param.upbnd")
    }
    if (is.null(lower))
	lower <- if(meth.has.bounds) asFinite(p.lowbnd + bound.eps*abs(p.lowbnd)) else -Inf
    if (is.null(upper))
	upper <- if(meth.has.bounds) asFinite(p.upbnd  - bound.eps*abs(p.upbnd )) else Inf

    logL <-
	if(traceOpt) {
	    LL <- function(param, u, copula) {
		r <- loglikCopula(param, u, copula)
		if(length(param) <= 1)
		    cat(sprintf("param=%14.9g => logL=%12.8g\n", param, r))
		else
		    cat(sprintf("param= %s => logL=%12.8g\n",
				paste(format(param, digits=9), collapse=", "), r))
		r
	    }
            if(need.finite) {
                nb <- length(body(LL))
                body(LL)[[nb]] <- do.call(substitute,
                                          list(body(LL)[[nb]],
                                               list(r = quote(asFinite(r)))))
            }
            LL
	} else if(need.finite) {
            function(param, u, copula) asFinite(loglikCopula(param, u, copula))
        } else
	    loglikCopula

    ## Maximize the likelihood
    fit <- optim(start, logL, lower=lower, upper=upper,
                 method = optim.method, control=control,
                 copula = copula, u = u)

    ## Check convergence of the fitting procedure
    freeParam(copula) <- fit$par
    loglik <- fit$val
    has.conv <- fit[["convergence"]] == 0
    if (is.na(estimate.variance))
        estimate.variance <- has.conv && ## for now:
	    (if(is(copula,"mixCopula"))
                 !is.null(wf <- attr(copula@w, "fixed")) && all(wf) else TRUE)
    if(!has.conv)
        warning("possible convergence problem: optim() gave code=",
                fit$convergence)

    ## Estimate the variance of the estimator
    var.est <- switch(method,
    "mpl" = {
        if(estimate.variance)
            var.mpl(copula, u) / nrow(u)
        else
            matrix(NA_real_, 0, 0) #matrix(NA_real_, q, q)
    },
    "ml" = {
        ## MM{FIXME}: This should be done only by 'vcov()' and summary() !
        if(estimate.variance) {
            fit.last <- optim(fit$par, # copula@@parameters
                              logL, lower=lower, upper=upper,
                              method=optim.method, copula=copula, u=u,
                              control=c(control, maxit=0), hessian=TRUE)
            vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
            if(is(vcov, "error")) {
                warning("Hessian matrix not invertible: ", vcov$message)
                matrix(NA_real_, 0, 0) #matrix(NA_real_, q, q)
            } else vcov ## ok
        } else matrix(NA_real_, 0, 0) #matrix(NA_real_, q, q)
    },
    stop("Wrong 'method'"))

    ## Return the fitted copula object
    new("fitCopula",
        estimate = fit$par,
        var.est = var.est,
        method = if(method=="mpl") "maximum pseudo-likelihood" else "maximum likelihood",
        loglik = loglik,
        fitting.stats = c(list(method=optim.method),
                          fit[c("convergence", "counts", "message")], control),
        nsample = nrow(u), call = call,
        copula = copula)
}


### Wrapper ####################################################################

##' @title Default fitCopula() Method -- *THE* user fitting function
##' @param copula The copula to be fitted
##' @param data The data in [0,1]^d for "mpl", "ml", "itau.mpl";
##'        for "itau", "irho", it can be in [0,1]^d or IR^d
##' @param method The estimation method
##' @param posDef A logical indicating whether pairwise estimated correlation
##'        matrices are turned into proper correlation matrices (via nearPD())
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param hideWarnings A logical indicating whether warnings from the
##'        underlying optimizations are suppressed
##' @param ... Additional arguments passed to auxiliary functions
##' @return The fitted copula object
fitCopula_dflt <- function(copula, data,
                           method = c("mpl", "ml", "itau", "irho", "itau.mpl"),
                           posDef = is(copula, "ellipCopula"),
                           start=NULL, lower=NULL, upper=NULL,
                           optim.method = optimMeth(copula, method, dim = d),
                           optim.control = list(maxit=1000),
                           estimate.variance = NA, hideWarnings = FALSE, ...)
{
    stopifnot(any(isFree(copula)), # GETR
	      is.list(optim.control) || is.null(optim.control))
    if(!is.matrix(data)) {
        warning("coercing 'data' to a matrix.")
        data <- as.matrix(data); stopifnot(is.matrix(data))
    }
    method <- match.arg(method)
    cl <- match.call()
    d <- ncol(data)
    if(method == "mpl" || method == "ml") { # "mpl" or "ml"
	if(is.function(optim.method)) ## << force(.) + flexibility for the user
	    optim.method <- optim.method(copula, method, dim=d)
        (if(hideWarnings) suppressWarnings else identity)(
        fitCopula.ml(copula, u=data, method=method,
                     start=start, lower=lower, upper=upper,
                     optim.method=optim.method, optim.control=optim.control,
                     estimate.variance=estimate.variance, call=cl, ...)
        )
    } else if(method == "itau" || method == "irho") { # "itau" or "irho"
        fitCopula.icor(copula, x=data, method=method,
                       estimate.variance=estimate.variance, call=cl, ...)
    } else { # "itau.mpl"
	if(!missing(optim.method))
	    warning(gettextf("method \"%s\" uses optimize(); hence '%s' is not used",
			     "itau.mpl", "optim.method"), domain=NA)
	if(!is.null(optim.control$trace)) {
	    warning(gettextf("method \"%s\" uses optimize(); hence '%s' is not used",
			     "itau.mpl", "optim.control = list(.., trace = *)"),
		    domain=NA)
	    warning("Rather use fitCopula(...., traceOpt = TRUE)")
	}
        (if(hideWarnings) suppressWarnings else identity)(
        fitCopula.itau.mpl(copula, u=data, posDef=posDef, lower=lower, upper=upper,
                           estimate.variance=estimate.variance, call=cl, ...) # <- may include 'tol' !
        )
    }
}

##' Default optim() method for  fitCopula() [and hence gofCopula(), xvCopula(),..]
optimMeth <- function(copula, method, dim) {
    ## Till Aug.2016, this was implicitly always "BFGS"
    cld <- getClass(class(copula))
    ellip <- extends(cld, "ellipCopula")
    if(ellip && copula@dispstr %in% c("ex", "ar1"))
        "L-BFGS-B"
    else if(extends(cld, "archmCopula")) {
	## if(nParam(copula, freeOnly=TRUE) == 1) # all 5 of our own families
	##     "Brent"
        ## "Brent" with "finite" bounds (almost Inf) --> quickly goes berserk
	## ## else if(class(copula) %in% archm.neg.tau && dim == 2)
	## ##     "Nelder-Mead" # for many theta, lik. L() = 0 (<==> log L() = -Inf)
	## else
        "L-BFGS-B"
    } else # "by default" -- was *the* only default till Aug. 2016:
        "BFGS"
    ## FIXME:  "Nelder-Mead" is sometimes better
}

setMethod( "fitCopula", signature("parCopula"), fitCopula_dflt)

## --> ./rotCopula.R :  "rotCopula"  has its own method
