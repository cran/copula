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

## FIXME: This has much in common with print.fitMvdc [ ./fitMvdc.R ]
## -----  For consistency, use common "helper" functions (instead of now: manually sync'ing)
print.fitCopula <- function(x, digits = max(3, getOption("digits") - 3), coefs = NULL,
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
	if(is.null(coefs)) coefs <- coef.fittedMV(x, SE = TRUE)
	printCoefmat(coefs, digits=digits, signif.stars=signif.stars,
		     na.print = "NA", ...)
    } else { # !showMore, default print() etc
	if(is.null(coefs)) coefs <- coef.fittedMV(x) # vector of "hat(theta)"
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

## NB: coef.fittedMV() does apply to <fitCopula>  --> ./Classes.R

## Keeping this back compatible... --> must return list(method=, loglik= ..):
summary.fitCopula <- function(object, orig=TRUE, ...) {
    structure(class = "summary.fitCopula",
              list(fitC = object,
                   method = object@method,
                   loglik = object@loglik,
                   convergence = object@fitting.stats[["convergence"]],
                   coefficients = coef.fittedMV(object, SE = TRUE, orig=orig)))
}

## Used both "as"  'print.summary.fitCopula()' and 'print.summary.fitMvdc()'
## Now print(summary(<fitCopula>)) *does* show a bit more than print(<fitCopula>)
printSummary.fittedMV <- function(x, ...) {
    print(x$fitC, coefs = x$coefficients, orig=x$orig, ..., showMore = TRUE)
    invisible(x)
}

setMethod("show", signature("fitCopula"),
	  function(object) print.fitCopula(object))

## Not really:  User should  use  dCopula(u, fitC @ copula, ..)  <<<<<<<<<
## setMethod("pCopula", signature("numeric", "fitCopula"),
##           function(u, copula, ...) pCopula(u, copula@copula, ...))
## setMethod("dCopula", signature("numeric", "fitCopula"),
##           function(u, copula, ...) dCopula(u, copula@copula, ...))
## setMethod("rCopula", signature("numeric", "fitCopula"),
##           function (n, copula, x, ...) rCopula(n, copula@copula, x, ...))

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

##' @title Computing the Log-Likelihood of the Given Copula   --> ../man/fitCopula.Rd
##' @param param *free* parameter (vector)
##' @param u The data in [0,1]^d
##' @param copula Copula object
##' @return Log-likelihood of the given copula at param given the data x
loglikCopula <- function(param = getTheta(copula), u, copula,
                         error = c("-Inf", "warn-Inf", "let-it-be"))
{
    if(!missing(param)) { # set the copula's parameter
        error <- match.arg(error)
        switch(error,
        "-Inf" = ,
        "warn-Inf" = {
            r <- tryCatch(freeParam(copula) <- param, error = function(e) e)
            if(inherits(r, "error")) {# param:  wrong-length || "out-of-bound" (TODO? look at error, maybe warn?)
                if(error == "warn-Inf") {
                    r$call <- sys.call()
                    warning(r)
                }
                return(structure(-Inf, reason = "param-setting error"))
            }
        },
        "let-it-be" =  # be fast in regular cases
            freeParam(copula) <- param,
        ## should never happen:
        stop("invalid 'error' argument (should not happen, please report!): ", error))
    }
    cop.param <- getTheta(copula, attr = TRUE) # freeOnly=TRUE;  attr=T : get bounds
    lower <- attr(cop.param, "param.lowbnd")
    upper <- attr(cop.param, "param.upbnd")
    admissible <- !any(is.na(cop.param) | cop.param > upper | cop.param < lower)
    if (admissible) {
        ## FIXME-JY: param range check is only a *part* of validity check
        ## MM: Requiring that copula's dCopula() method *must*, return() a number(also -Inf, NA) here:
        sum(dCopula(u, copula=copula, log=TRUE, checkPar=FALSE))
    } else structure(-Inf, reason = "inadmissible param")
}

##' not exported; used in  fitCopStart() and getIniParam()'s default method:
getIni_itau <- function(copula, u, default, classDef = getClass(class(copula)), ...)
{
    .par.df <- has.par.df(copula, classDef)
    start <- fitCopula.icor(if(.par.df) as.df.fixed(copula, classDef = classDef) else copula,
                            x=u, method="itau", estimate.variance=FALSE,
                            warn.df=FALSE, ...)@estimate # fitCopula.icor(, method="itau")
    if(.par.df) # add starting value for 'df'
        start <- c(start, getdf(copula))
    ll <- loglikCopula(start, u=u, copula=copula)
    if(is.finite(ll))
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
        else if(!is.na(ll) && ll == +Inf)
            start # e.g. for perfectly correlated data
        else
            default
    }
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
    if(!is.null(FUN <- selectMethod("getIniParam", signature=clc, optional = TRUE)))
        FUN(copula, u, default, named=TRUE, ...)
    else if(hasMethod("iTau", clc))
        getIni_itau(copula, u, default, classDef=getClass(clc), ...)
    else
        default
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
    dim <- dim(copula)
    dd <- dim * (dim - 1) / 2
    free <- isFree(copula)
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
    dim <- dim(copula)
    dd <- dim * (dim - 1) / 2
    xmat <-
    if (!is(copula, "ellipCopula")) # one-parameter non-elliptical copula
	if((n.th <- nParam(copula, freeOnly=TRUE)) == 1L)
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
    free <- isFree(copula) ## the last one is for df and not needed
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
    dim <- dim(cop)
    dd <- dim * (dim - 1) / 2
    free <- isFree(cop)
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
	    return(matrix(numeric(), 0, 0)) # instead of potentiall large (n x dd)
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
    if (!isEll && !hasMethod(dlogcdu, clc)) {
	warning(msg); return(ans)
    }

    ## If df.fixed = FALSE, Jscore() cannot be computed
    if(has.par.df(cop, ccl, isEll)) {
        cop <- as.df.fixed(cop, classDef = ccl)
        warning(
     "the covariance matrix of the parameter estimates is computed as if 'df.fixed = TRUE' with df = ",
                getdf(cop))
        ans[-p, -p] <- var(t(Jscore(cop, u, method = "mpl")))
        ans
    }
    else if(is(cop, "mixCopula") && any(isFreeP(cop@w))) {
                                        # mixCopula with estimated weights
        ans[-p, -p] <- var(t(Jscore(cop, u, method = "mpl")))
        ans[p, ] <- ans[,p] <- 0.
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
    stopifnot(any(free <- isFree(copula)))
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
    } else matrix(numeric(), 0, 0)
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
    stopifnot(any(free <- isFree(copula)))
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    ## Note: We require dispstr = "un" as it is a) the method of Mashal, Zeevi (2002) and
    ##       b) otherwise would require the use of .lm.fit which can crash for large d!
    ## MM: The above seems wrong conceptually, "un" being ill-posed for large d
    if(!(is(copula, "tCopula") && copula@dispstr == "un"))
        stop("method \"itau.mpl\" is only applicable for \"tCopula\" with 'dispstr=\"un\"'")
    if(copula@df.fixed) stop("Use method=\"itau\" for 'tCopula' with 'df.fixed=TRUE'")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (dim(copula) != d)
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
            numTrace <- (is.numeric(traceOpt) && (traceOpt <- as.integer(traceOpt)) > 0)
            ## if(numTrace) "print every traceOpt iteration"
	    function(Inu) {
		r <- loglikCopula(c(P, df=1/Inu)[free], u=u, copula=copula) # IK: [free]
                if(numTrace) iTr <<- iTr+1L
                if(!numTrace || iTr == 1L || (iTr %% traceOpt == 0L))
                    cat(sprintf("%s1/nu=%14.9g => nu=%9.4f; logL=%12.8g\n",
                                if(numTrace) sprintf("%3d: ", iTr) else "",
                                Inu, 1/Inu, r))
		r
	    }
	} else
	    function(Inu) loglikCopula(c(P, df=1/Inu)[free], u=u, copula=copula) # IK: [free]

    if(traceOpt && numTrace) iTr <- 0L
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
    var.est <- matrix(numeric(), 0, 0) # length 0  <==> not-estimated / not-available

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

##' @title Maximum Likelihood Estimator for Copulas -- not exported; called from fitCopula()
##' @param copula The copula to be fitted
##' @param u The data in [0,1]^d (for method="ml", this needs to be true copula data;
##'        for method="mpl", this can be parametrically or non-parametrically estimated
##'        pseudo-observations)
##' @param method string, either "mpl" (default) or "ml"
##' @param start The initial value for optim()
##' @param lower The vector of lower bounds for optim()
##' @param upper The vector of upper bounds for optim()
##' @param optim.method The optimization method for optim()
##' @param optim.control optim()'s control parameter
##' @param estimate.variance A logical indicating whether the estimator's
##'        variance shall be computed
##' @param bound.eps A small quantity denoting an eps for staying away from
##'        the theoretical parameter bounds
##' @param call the call, typically passed from caller
##' @param traceOpt logical or positive integer indicating if the object function should be traced during minimization
##' @param need.finite logical indicating if the optimizer needs \emph{finite} obj.fn values
##' @param finiteLARGE large positive number to replace \code{Inf} when \code{need.finite}
##'                    is true or the \code{optim.method} is that needs box constraint bounds.
##' @param ... Additional arguments (currently with no effect)
##' @return The fitted copula object
fitCopula.ml <- function(copula, u, method=c("mpl", "ml"), start, lower, upper,
                         optim.method, optim.control, estimate.variance,
                         bound.eps=.Machine$double.eps^0.5, call, traceOpt = FALSE,
			 need.finite = any(optim.method == c("Brent", "L-BFGS-B", "BFGS")),
			 finiteLARGE = .Machine$double.xmax,
			 ...)
{
    q <- nParam(copula, freeOnly=TRUE)
    stopifnot(q > 0L,
              is.list(optim.control) || is.null(optim.control),
              is.character(optim.method), length(optim.method) == 1,
              is.finite(finiteLARGE), length(finiteLARGE) == 1, finiteLARGE > 0)
    chk.s(...) # 'check dots'
    if(any(u < 0) || any(u > 1))
        stop("'u' must be in [0,1] -- probably rather use pobs(.)")
    stopifnot(is.numeric(d <- ncol(u)), d >= 2)
    if (dim(copula) != d)
        stop("The dimension of the data and copula do not match")
    if(missing(call)) call <- match.call()
    if(is.null(start))
        start <- fitCopStart(copula, u=u)
    if(anyNA(start)) stop("'start' contains NA values")

    if(ismixC <- is(copula, "mixCopula")) {
        is.freeW <- isFreeP(copula@w)
        nfree.w <- sum(is.freeW) # number of free weights
        if(nfree.w == 1) {
            warning("only one mixture weight is free; this is illogical, setting all to fixed")
            attr(copula@w, "fixed") <- TRUE
            nfree.w <- sum(is.freeW <- isFreeP(copula@w)) ## = 0; also update 'is.freeW'
        }
        ## ------ (keep these checks for now, but they never triggered  ---------
        i.free.w <- which(is.freeW)
        if(nfree.w != length(i.free.w))
            stop(gettextf("fitCopula.ml(<mixCopula>): nfree.w = %d != length(i.free.w) = %d",
                          nfree.w, length(i.free.w)), domain=NA)
        ## number of free copula parameters
        ## nfreecop <- sum(c(unlist(lapply(copula@cops, isFree))))
        ## nfreecop <- sum(vapply(copula@cops, function(C) nFree(C@parameters), 0L))
        ## nfreecop <- sum(vapply(copula@cops, nFree, 0L))
        ## if(nfreecop+nfree.w != q) ## no longer true
        ##     stop(gettextf("fitCopula.ml(<mixCopula>): (nfreecop+nfree.w) = %d != q = %d",
        ##                   nfreecop+nfree.w, q), domain=NA)
        ## ------
        ## when there are no free weights, don't need transformations:
        ismixC <- nfree.w > 0 # and it is  >= 2
    }   ##====

    if(q != length(start)) ## in the "simple" (non-mixture) case:
        stop(gettextf(
	"The lengths of 'start' (= %d) and the number of free copula parameters (=%d) differ",
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
        ## _Correct_bounds_  >>>> not-yet__see_*more*_problems__gofCopula(normalCopula()..)
        ## if (is.null(lower)) {
        ##     notI <- !is.infinite(p.lowbnd)
        ##     p.lowbnd[notI] <- p.lowbnd[notI] + bound.eps*abs(p.lowbnd[notI])
        ## }
        ## if (is.null(upper)) {
        ##     notI <- !is.infinite(p.upbnd)
        ##     p.upbnd [notI] <- p.upbnd [notI] - bound.eps*abs(p.upbnd [notI])
        ## }
    }
    if (is.null(lower))
        ## _Correct_bounds_ lower <- if(meth.has.bounds) asFinite(p.lowbnd) else -Inf
        lower <- if(meth.has.bounds) asFinite(p.lowbnd + bound.eps*abs(p.lowbnd)) else -Inf

    if (is.null(upper))
        ## _Correct_bounds_ upper <- if(meth.has.bounds) asFinite(p.upbnd ) else Inf
        upper <- if(meth.has.bounds) asFinite(p.upbnd  - bound.eps*abs(p.upbnd )) else Inf

    if(ismixC) {
        ## m <- length(is.freeW) # == length(copula@w)
        sumfreew <- 1 - sum(copula@w[!is.freeW]) ## sum of free weights

        ## NB: start *now* contains m weights in the original m-simplex (if all 'w' are free),
        ##     respectively  nfree.w  weights in the nfree.w-simplex (with sum 'sumfreew')
        ##     if *some* w_j are fixed (but at least two are free!):
        ## Go to the (m-1) - dimensional transformed ("lambda") space for logL() & optim(.),
        ## see also loglikMixCop() :
        ## q == length(start)
        l  <- q - 1L
        i1 <- q - (nfree.w - 1L)
        in1 <- i1:l
        in2 <- i1:q
        start[in1] <- clr1(start[in2] / sumfreew)
        length(start) <- l # cutoff last one
        ## similarly transform lower[], upper[] for the 'w' -> 'l' transformation:
        if(meth.has.bounds) {
            ## --- fails --- as "weights don't add to 1" ---
            ## lower[in1] <- clr1(lower[in2] / sumfreew); length(lower) <- l
            ## upper[in1] <- clr1(upper[in2] / sumfreew); length(upper) <- l
            lower[in1] <- -Inf; length(lower) <- l
            upper[in1] <- +Inf; length(upper) <- l
        }
        ##
        rm(l, i1, q) # not used below, whereas (in1, in2, sumfreew) *are* needed

        ## the last nfree.w   elements of the vector 'start' correspond to the
        ##          nfree.w-1 transformed weights '\lambda's
        ## transform them back to the nfree.w 'w's:
        loglikMixCop <- function(param, u, copula) { # + (sumfreew, in1, in2)  "global"s
            ## extend 'param' length by 1 with weights of length 'nfree.w' :
            param[in2] <- clr1inv(param[in1]) * sumfreew
            loglikCopula(param, u, copula)
        }
    } # if(ismixC)

    logL <-
        if(traceOpt) {
            numTrace <- (is.numeric(traceOpt) && (traceOpt <- as.integer(traceOpt)) > 0)
            ## if(numTrace) "print every traceOpt iteration"
            LL <- function(param, u, copula) {
                r <- if(ismixC)
                          loglikMixCop(param, u, copula)
                     else loglikCopula(param, u, copula)
                if(numTrace) iTr <<- iTr+1L
                if(!numTrace || iTr == 1L || (iTr %% traceOpt == 0L)) {
                    chTr <- if(numTrace) sprintf("%3d:", iTr) else ""
                    if(length(param) <= 1)
                        cat(sprintf("%s param=%14.9g => logL=%12.8g\n", chTr, param, r))
                    else
                        cat(sprintf("%s param= %s => logL=%12.8g\n", chTr,
                                    paste(format(param, digits=9), collapse=", "), r))
                }
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
            if(ismixC)
                 function(param, u, copula) asFinite(loglikMixCop(param, u, copula))
            else function(param, u, copula) asFinite(loglikCopula(param, u, copula))
        } else if(ismixC) {
            ## loglikCopula with inverse transformed lambdas, i.e., the original 'w's:
            loglikMixCop
        } else
            loglikCopula

    if(traceOpt && numTrace) iTr <- 0L
    ## Maximize the (log) likelihood
    fit <- optim(start, logL, lower=lower, upper=upper,
                 method = optim.method, control=control,
                 copula = copula, u = u)

    fitpar <- fitp0 <- fit$par
    if(ismixC) { ## Transform the optimized '\lambda's back
        ## defined (in1, in2)  above
        fitpar[in2] <- clr1inv(fitp0[in1]) * sumfreew
    }

    freeParam(copula) <- fitpar
    ## Check convergence of the fitting procedure
    loglik <- fit$val
    has.conv <- fit[["convergence"]] == 0
    if (is.na(estimate.variance))
        estimate.variance <- has.conv && ## for now:
            (if(ismixC)
                 !is.null(wf <- attr(copula@w, "fixed")) && all(wf) else TRUE)
    if(!has.conv)
        warning("possible convergence problem: optim() gave code=",
                fit$convergence)

    ## Estimate the variance of the estimator
    Var0  <- matrix(numeric(), 0, 0)
    VarNA <- matrix(NA_real_, length(start), length(start))
    var.est <- switch(
        method,
        "mpl" = {
            if(estimate.variance) ## FIXME?  use tryCatch() as in "ml" case ??
                var.mpl(copula, u) / nrow(u)
            else
                Var0
        },
        "ml" = {
            ## FIXME: should be done additionally / only by 'vcov()' and summary()
            if(estimate.variance) {
                if(traceOpt)
                    cat("--- estimate.variance=TRUE ==> optim(*, hessian=TRUE) iterations ---\n")
                fit.last <- tryCatch(
                    optim(fitp0, # typically == freeParam(copula)
                          logL, lower=lower, upper=upper,
                          method=optim.method, copula=copula, u=u,
                          control=c(control, maxit=0), hessian=TRUE),
                    error = function(e) e)
                if(inherits(fit.last, "error")) {
                    msg <- .makeMessage("optim(*, hessian=TRUE) failed: ", fit.last$message)
                    warning(msg)
                    structure(VarNA, msg = msg)
                } else {
                    vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
                    if(inherits(vcov, "error")) {
                        msg <- .makeMessage("Hessian matrix not invertible: ", vcov$message)
                        warning(msg)
                        structure(VarNA, msg = msg)
                    } else vcov ## ok
                }
            } else Var0
        },
        stop("Wrong 'method'"))
    ## Return the fitted copula object
    new("fitCopula",
        estimate = fitpar,
        var.est = var.est,
        method = if(method=="mpl") "maximum pseudo-likelihood" else "maximum likelihood",
        loglik = loglik,
        fitting.stats = c(list(method=optim.method),
                          fit[c("convergence", "counts", "message")],
                          control,
                          if(ismixC)
                              list(paramTrafo =
                                       list(kind = "mixCop-clr1",
                                            sumfreew=sumfreew, in1=in1, in2=in2))),
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
    stopifnot(any(isFree(copula)),
	      is.list(optim.control) || is.null(optim.control))
    if(!is.matrix(data)) {
        warning("coercing 'data' to a matrix.")
        data <- as.matrix(data); stopifnot(is.matrix(data))
    }
    method <- match.arg(method)
    cl <- match.call()
    if(cl[[1]] == quote(.local)) ## fix it up -- but how?  This fails:
        {} ## FIXME: cl <- match.call(sys.function(sys.parent()))
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
    else if(extends(cld, "archmCopula") || extends(cld, "mixCopula") ) {
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
