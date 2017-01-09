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

setClass("fitMvdc", contains = "fittedMV", #-> ./Classes.R
         slots = c(mvdc = "mvdc")
	 ## FIXME , validity = function(object) TRUE
	 )

paramNmsMvdc <- function(mv) c(margpnames(mv), paramNames(mv@copula))
## not exported, but used ..
setMethod("paramNames", "fitMvdc", function(x) paramNmsMvdc(x@mvdc))

## FIXME: This has much in commong with print.fitCopula [ ./fitCopula.R ]
## -----  For consistency, use common "helper" functions (instead of now: manually sync'ing)
print.fitMvdc <- function(x, digits = max(3, getOption("digits") - 3),
                          signif.stars = getOption("show.signif.stars"), ...,
                          showMore = FALSE)
{
    cat("Call: ", formatCall(x@call, class(x)), "\n", sep = "")
    cop <- x@mvdc@copula
    d <- dim(cop) # not slot; e.g. for rotCopula
    cat(sprintf(
	"Maximum Likelihood estimation based on %d %d-dimensional observations.\n",
	x@nsample, d))
    ## FIXME show more via printCopula() utility; but do *not* show the parameters
    cat("Copula: ", class(cop), "\n")
    if(showMore) { ## coefficient *matrix*, incl. std.errs
        coefs <- coef.fittedMV(x, SE = TRUE)
        printCf <- function(ind)
	    printCoefmat(coefs[ ind, , drop=FALSE],
			 digits = digits, signif.stars = signif.stars,
			 signif.legend=FALSE, ## in any case
			 na.print = "NA", ...)
    } else { # !showMore, default print() etc
        coefs <- coef.fittedMV(x) # vector of "hat(theta)"
        printCf <- function(ind) print(coefs[ind], digits=digits, ...)
    }
    marNpar <- lengths(x@mvdc@paramMargins)# or  vapply(x@mvdc@paramMargins, nFree, 1L)
    idx2 <- cumsum(marNpar)
    idx1 <- idx2 - marNpar + 1
    margid <- x@mvdc@marginsIdentical
    if (sum(marNpar) > 0) { ## sometimes there is no marginal params
	if(margid){
	    cat("Identical margins:\n")
            printCf(idx1[1]:idx2[1])
	}
	else for (i in 1:d) {
                 cat("Margin", i, ":\n")
                 printCf(idx1[i]:idx2[i])
             }
    }
    cat(if(showMore) describeCop(cop, "short") # as 'parameters' are separate
        else "Copula:", "\n")
## or	else paste("Copula:", class(x)), "\n")
    copParIdx <- seq_len(nFree(cop@parameters))
    printCf(copParIdx + if(margid) marNpar[1] else sum(marNpar))
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
    invisible(x)
}

summary.fitMvdc <- function(object, ...) {
  structure(class = "summary.fitMvdc",
	    list(fitC = object,
                 method = object@method,
                 loglik = object@loglik,
		 convergence = object@fitting.stats[["convergence"]],
                 coefficients = coef.fittedMV(object, SE = TRUE)))
}

## NB: print.summary.fitMvdc() via  printSummary.fittedMV() [ ./fitCopula.R ]

setMethod("show", signature("fitMvdc"), function(object) print.fitMvdc(object))


################################################################################

##' @title Set / Update the free (non-fixed) parameters of an MVDC
##' @param mvdc an \code{mvdc} object
##' @param param numeric vector of free parameters
##' @return modified \code{mvdc} object
setMvdcPar <- function(mvdc, param, noCheck = FALSE) {
    marNpar <- lengths(mvdc@paramMargins)# or  vapply(mvdc@paramMargins, nFree, 1L)
    idx2 <- cumsum(marNpar)
    idx1 <- idx2 - marNpar + 1
    margid <- mvdc@marginsIdentical
    d <- dim(mvdc@copula)

    for (i in 1:d) {
        if (marNpar[i] > 0) {
            ## parnames <- mvdc@paramMargins[[i]]
            k <- if(margid) 1 else i
            par <- param[idx1[k]: idx2[k]]
            ## names(par) <- parnames
            ## mvdc@paramMargins[i] <- as.list(par)
            for (j in 1:marNpar[i]) mvdc@paramMargins[[i]][j] <- par[j]
        }
    }
    mvdc@copula <-
        setTheta(mvdc@copula,
                 value =
                     if (idx2[d] == 0) param # no marginal parameters
                     else param[- (if(margid) 1:idx2[1] else 1:rev(idx2)[1])],
                 noCheck = TRUE)
    mvdc
}

loglikMvdc <- function(param, x, mvdc) {
  mvdc <- setMvdcPar(mvdc, param, noCheck=TRUE)
  tryCatch(
      ## log likelihood :
      sum(log(dMvdc(x, mvdc))),
      error = function(e) {
	  warning("error in loglik computation: ", conditionMessage(e))
	  (-Inf) # was NaN
      })
}

fitMvdc <- function(data, mvdc, start,
                    optim.control=list(), method="BFGS",
                    lower = -Inf, upper = Inf,
                    estimate.variance = fit$convergence == 0, hideWarnings=TRUE)
{
    copula <- mvdc@copula
    if (copula@dimension != ncol(data))
        stop("The dimensions of the data and copula do not match.")
    cl <- match.call()
    marNpar <- lengths(mvdc@paramMargins)# or  vapply(mvdc@paramMargins, nFree, 1L)
    margid <- mvdc@marginsIdentical
    q <- length(start)
    if(q != nParam(copula, freeOnly=TRUE) + (if(margid) marNpar[1] else sum(marNpar)))
	stop("The lengths of 'start' and mvdc parameters do not match.")
    mvdCheckM(mvdc@margins, "p")
    control <- c(optim.control, fnscale=-1)
    control <- control[ !vapply(control, is.null, NA)]

    ## messageOut may be used for debugging
    (if(hideWarnings) suppressWarnings else identity)(
    fit <- optim(start, loglikMvdc,
		 ## loglikMvdc args:
		 mvdc=mvdc, x=data,
		 ## optim args:
		 method=method, control=control, lower=lower, upper=upper)
    )

    if (fit$convergence > 0)
	warning("possible convergence problem: optim gave code=", fit$convergence)
    loglik <- fit$val
    param <- fit$par

    varNA <- matrix(NA_real_, q, q)
    var.est <- if (estimate.variance) {
	fit.last <- optim(param, loglikMvdc,
			  ## loglikMvdc args :
			  mvdc=mvdc, x=data,
			  ## optim args:
			  method = method, ## one final step, computing Hessian :
			  control=c(control, maxit = 1), hessian=TRUE)
        fit$counts <- fit$counts + fit.last$counts
	vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
	if(is(vcov, "error")) {
	    warning("Hessian matrix not invertible: ", vcov$message)
	    varNA
	} else vcov ## ok
    } else varNA

    new("fitMvdc",
	estimate = param,
	var.est = var.est,
	loglik = loglik,
	method = method,
	fitting.stats = c(fit[c("convergence", "counts", "message")], control),
	nsample = nrow(data),
	call = cl,
	## this contains 'copula':
	mvdc = setMvdcPar(mvdc, param))
}
