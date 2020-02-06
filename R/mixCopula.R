##' A list of parCopula's
setClass("parClist",  contains = "list",
	 validity = function(object) {
	     if(!length(object))
		 return("empty parCopula lists are not valid")
	     if(!all(lengths(cls <- lapply(object, class)) == 1L))
		 return("parCopula list must have parCopula components")
	     is.pC <- vapply(unlist(cls), extends, NA, class2 = "parCopula")
	     if(!all(is.pC))
		 return(paste("components", paste(which(!is.PC), collapse=", "),
			      "do not inherit from \"parCopula\""))
	     dims <- vapply(object, dim, 1)
	     if(any(ne.d <- dims[1] != dims))
		 return(paste("copula dimensions differ", dims[1], "!=",
			      dims[match(FALSE, ne.d)]))
	     ## else
	     TRUE
	 })

##' Mixture Weights (not yet exported)
setClass("mixWeights", contains = "numeric",
	 validity = function(object) {
	     if(any(object < 0))
		 return("mixture weights must not be negative")
	     s <- sum(object)
	     if(abs(s - 1) > 64*.Machine$double.eps)
		 return(paste("mixture weights do not add to one, but to 1 -",
			      formatC(1-s)))
	     TRUE
	 })

setAs("numeric", "mixWeights", function(from) {
    if(any(neg <- from < 0)) {
        warning("negative mixture weights set to 0")
        from[neg] <- 0
    }
    new("mixWeights", from / sum(from)) # coerce sum(.) == 1
})


##' A Mixture of Copulas
setClass("mixCopula", contains = "parCopula",
	 slots = c("w" = "mixWeights", "cops" = "parClist"),
	 validity = function(object) {
	     if((n1 <- length(object@w)) != (n2 <- length(object@cops)))
		 return(paste("length of weights and copula list differ:",
			      paste(n1, "!=", n2)))
	     TRUE
	 })

##' A Mixture of Explicit Copulas
setClass("mixExplicitCopula", contains = "mixCopula",
         slot = c("exprdist" = "expression")
         )


##################################################################################
### Constructor of mixed copulas
##################################################################################

##' Creates a mixCopula object or a mixExplicitcopula object
##'
##' @title Creates a mixcopula or mixExplicitcopula object
##' @param coplist a list of copulas of the same dimension
##' @param w a numeric vector of the same length as coplist for the mixing weights
##' @return a new mixCopula or mixExplicitcopula object
##'
mixCopula <- function(coplist, w = NULL) {
    stopifnot(is.list(coplist))
    if((m <- length(coplist)) == 0)
	stop("a mixCopula must have at least one copula component")
    if (is.null(w)) # default: equal weights
	w <- rep.int(1/m, m)
    allExplicit <- all(vapply(coplist, isExplicit, NA))

    ## now the validity methods check:
    if (!allExplicit)
        new("mixCopula",
            w = as(w, "mixWeights"),
            cops = as(coplist, "parClist"))
    else { ## mixExplicitCopula
        mixCdf <- mixPdf <- parse(text = paste0("pcdf", 1L:m, " * ", "w", 1L:m,
                                                collapse = " + "))
        exprdist <- c(cdf = mixCdf, pdf = mixPdf)
        cdfL <- lapply(coplist, function(x) x@exprdist$cdf)
        pdfL <- lapply(coplist, function(x) x@exprdist$pdf)
        ListExpr <- function(nm1, nm2) # ((uglyness in one place))
            parse(text = paste0("list(", paste0(nm1, "= quote(",nm2,")", collapse= ", "), ")"))
        for (i in 1:m) {
            ## original 6 basic families have alpha in expressions
            ## cdfL[[i]] <- do.call(substitute, list(cdfL[[i]], list(alpha = quote(param))))
            ## pdfL[[i]] <- do.call(substitute, list(pdfL[[i]], list(alpha = quote(param))))
            ## rename the parameters with a prefix of 'm<i>.'
            oldParNames <- names(getTheta(coplist[[i]], freeOnly=FALSE, named=TRUE))
            npar <- length(oldParNames)
            if (npar > 0) {
                prefix <- paste0("m", i, ".")
                newParNames <- paste0(prefix, oldParNames)
                rep.l <- ListExpr(oldParNames, newParNames)
                cdfL[[i]] <- do.call(substitute, list(cdfL[[i]], eval(rep.l)))
                pdfL[[i]] <- do.call(substitute, list(pdfL[[i]], eval(rep.l)))
            }
        }
        cdfL <- as.expression(cdfL)
        pdfL <- as.expression(pdfL)
        pcdfs <- paste0("pcdf", 1:m)
        cdf.repl <- ListExpr(pcdfs, cdfL)
        pdf.repl <- ListExpr(pcdfs, pdfL)
        ## why this does not work? what happened when they were put together with c?
        ## mixCdf <- do.call(substitute, list(mixCdf, eval(cdf.repl)))
        ## mixPdf <- do.call(substitute, list(mixPdf, eval(pdf.repl)))
        mixCdf <- do.call(substitute, list(exprdist$cdf, eval(cdf.repl)))
        mixPdf <- do.call(substitute, list(exprdist$pdf, eval(pdf.repl)))

        cdf <- as.expression(mixCdf)
        cdf.algr <- deriv(cdf, "nothing")
        pdf <- as.expression(mixPdf)
        pdf.algr <- deriv(pdf, "nothing")
        exprdist <- c(cdf = cdf, pdf = pdf)
        attr(exprdist, "cdfalgr") <- cdf.algr
        attr(exprdist, "pdfalgr") <- pdf.algr

        new("mixExplicitCopula",
            w = as(w, "mixWeights"),
            cops = as(coplist, "parClist"),
            exprdist = exprdist)
    }
}

##################################################################################
### Basic methods
##################################################################################

## dimension
setMethod("dim", signature("mixCopula"), function(x) dim(x@cops[[1]]))

## logical indicating which parameters are free
setMethod("isFree", signature("mixCopula"), function(copula)
    c(unlist(lapply(copula@cops, isFree)), isFreeP(copula@w))) # FIXME reparametrize 'w'

## number of (free / all) parameters :
setMethod("nParam", signature("mixCopula"), function(copula, freeOnly=FALSE)
    sum(vapply(copula@cops, nParam, 1, freeOnly=freeOnly)) +
    (if(freeOnly) nFree else length)(copula@w))
## FIXME reparametrize 'w'
## pmax(0, (if(freeOnly) nFree else length)(copula@w) - 1L))
##                                                     ==== "first" = 1 - sum(<others>)

## parameter names for freeOnly parameters
setMethod("paramNames", signature("mixCopula"), function(x) {
    cops <- x@cops
    ## m <- length(cops)
    nj <- vapply(cops, nParam, 1, freeOnly=TRUE)
    ic <- seq_along(cops)# 1:m
    c(unlist(lapply(ic,
             function(i) if(nj[i] > 0L)
                             paste0("m", i, ".", paramNames(cops[[i]])))),
      paste0("w", ic)[isFreeP(x@w)]) # FIXME reparametrize: -> theta
})

## get parameters -- used in loglikCopula() [in ./fitCopula.R ]
setMethod("getTheta", "mixCopula",
function(copula, freeOnly = TRUE, attr = FALSE, named = attr) {
    parC <- lapply(copula@cops, getTheta,
                   freeOnly=freeOnly, attr=attr, named=named)
    w <- copula@w
    if(freeOnly) {
	wFree <- isFreeP(w)
	w <- w[wFree]
    }
    if(named) {
        ic <- seq_along(parC) ## not w because w can be free only
        for(i in ic[lengths(parC) > 0])
            names(parC[[i]]) <- paste0("m", i, ".", names(parC[[i]]))
        attr(w, "names") <- paste0("w", ic)[if(freeOnly) wFree else TRUE]
    }
    ## FIXME re-parametrize 'w' a la nor1mix:: (??)
    if(attr) { # more structured result
        ## need attr(*, "param.(up|low)bnd") for for loglikCopula()
        lowb <- lapply(parC, attr, "param.lowbnd")
        uppb <- lapply(parC, attr, "param.upbnd")
        structure(c(unlist(parC),                    w),
                  param.lowbnd = c(unlist(lowb), rep(0, length(w))),
                  param.upbnd  = c(unlist(uppb), rep(1, length(w))))
    } else {
        c(unlist(parC), w) # FIXME reparametrize: -> theta
    }
})

## set free parameters [hidden, must be fast for loglikCopula <= fitCopula.ml]:
setMethod("freeParam<-", signature("mixCopula", "numeric"),
function(copula, value) {
    cops <- copula@cops
    m <- length(cops)
    nj <- vapply(cops, nParam, 1, freeOnly=TRUE)
    ## FIXME re-parametrize a la nor1mix::nM2par / .par2nM
    ## ----- i.e. value would only contain  qlogis(w[-1])  !!
    nw <- sum(iF.w <- isFreeP(w <- copula@w))
    if (sum(nj) + nw != length(value))
        stop(sprintf(
            "length(value) = %d  !=  %d, the number of free parameters",
            length(value), sum(nj) + nw))
    n. <- 0L
    for(j in seq_len(m)) # set parameters for component copulas
        if ((n.j <- nj[j])) {
            freeParam(cops[[j]]) <- value[n.+ seq_len(n.j)]
            n. <- n. + n.j
        }
    if(n.)
        copula@cops <- cops
    if(nw) { ## Ensuring that w sums to one (!)
	I. <- 1 - sum(w[!iF.w]) # == target sum(<free w>)
	w. <- value[n. + seq_len(nw)]
	if(any(w. < 0)) stop("mixCopula weights must not become negative")
	copula@w[iF.w] <- I. * w. / sum(w.)
    }
    if(getOption("copula:checkMix", FALSE))
	stopifnot(validObject(copula))
    copula
})

## set parameters
setMethod("setTheta", "mixCopula",
          function(x, value, na.ok = TRUE, noCheck = FALSE, freeOnly = TRUE,
                   treat.negative = c("set.0", "warn.set0", "stop"), ...) {
    stopifnot(is.numeric(value) | (ina <- is.na(value)))
    if(any(ina)) {
        if(!na.ok) stop("NA value, but 'na.ok' is not TRUE")
        ## vectorized (and partial) value <- NA_real_
        if(!is.double(value)) storage.mode(value) <- "double"
    }
    cops <- x@cops
    m <- length(cops)
    nj <- vapply(cops, nParam, 1, freeOnly=freeOnly)
    ## FIXME re-parametrize a la nor1mix::nM2par / .par2nM
    ## ----- i.e. value would only contain  qlogis(w[-1])  !!
    nw <- length(iF.w <- isFreeP(w <- x@w))
    if (sum(nj) + nw != length(value))
        stop(sprintf(
            "length(value) = %d  !=  %d, the number of free parameters",
            length(value), sum(nj) + nw))
    n. <- 0L
    for(j in seq_len(m))
        if ((n.j <- nj[j])) {
            freeParam(cops[[j]]) <- value[n.+ seq_len(n.j)]
            n. <- n. + n.j
        }
    if(n.)
        x@cops <- cops
    if(nw) { ## Ensuring that w sums to one
	I. <- 1 - sum(w[!iF.w]) # == target sum(<free w>)
	w. <- value[n. + seq_len(nw)]
	if(any(wN <- w. < 0)) {
	    switch(match.arg(treat.negative),
		   "stop" =
		       stop("mixCopula weights must not become negative"),
		   "warn.set0" =
		       warning("negative mixCopula weights set to 0"))
	    ## if not stopped above,
	    w.[wN] <- 0
	    ## if(sum(w.) == 0) ... what exactly?
        }
	x@w[iF.w] <- I. * w. / sum(w.)
    }
    if(!noCheck)
        stopifnot(validObject(x))
    x
})

## set or modify "fixedness" of parameters
setMethod("fixedParam<-", signature("mixCopula", "logical"),
function(copula, value) {
    stopifnot(length(value) %in% c(1L, nParam(copula)))
    if (anyNA(getTheta(copula, freeOnly = FALSE)[value])) stop("Fixed parameters cannot be NA.")
    cops <- copula@cops
    ic <- seq_along(cops) # "1:m"
    nj <- vapply(cops, nParam, 1, freeOnly=FALSE)
    if (identical(value, FALSE) || !any(value)) {
        for(j in ic)
            if (nj[j] > 0L) fixedParam(copula@cops[[j]]) <- FALSE
        attr(copula@w, "fixed") <- NULL
    } else if (identical(value, TRUE) || all(value)) {
        for(j in ic)
            if (nj[j] > 0L) fixedParam(copula@cops[[j]]) <- TRUE
        attr(copula@w, "fixed") <- TRUE
    } else { # "typically", some fixed, some free:
        n. <- 0L
        for (j in ic) {
            n.j <- nj[j]
            if (n.j > 0L) fixedParam(copula@cops[[j]]) <- value[n. + seq_len(n.j)]
            n. <- n. + n.j
        }
        attr(copula@w, "fixed") <-
            if (!any(vw <- value[n. + ic])) NULL
            else if (all(vw)) TRUE else vw
    }
    copula
})


## Describe copula
setMethod(describeCop, c("mixCopula", "character"), function(x, kind, prefix="", ...) {
    m <- length(x@w)
    c1 <- paste0(prefix, "mixCopula from ", m, " components")
    if(kind == "very short")
        return(c1)
    ## else
    dC <- vapply(x@cops, describeCop, "", kind="very short", prefix=" ")
    paste0(c1, "\n",
	   dputNamed(if(kind == "short") unname(vapply(dC, abbreviate, "")) else dC),
	   "  with weights:\n", dputNamed(x@w))
})

## The copula
pMixCopula <- function(u, copula, ...) {
    as.vector(
	vapply(copula@cops, pCopula, FUN.VALUE=numeric(nrow(u)), u=u, ...)
	%*%
	copula@w)
}

setMethod("pCopula", signature("matrix",  "mixCopula"), pMixCopula)

## The density
dMixCopula <- function(u, copula, log = FALSE, ...) {
    n <- nrow(u) ## stopifnot(is.numeric(n))
    fu <- vapply(copula@cops, dCopula, FUN.VALUE=numeric(n), u=u, log=log, ...)
    if(n == 1L) dim(fu) <- c(n, length(fu))# vapply() gave p-vector instead of (1,p) matrix
    w <- as.numeric(copula@w)
    if(log)
	lsum(t(fu) + log(w)) ## = log(colSums(exp(t(fu) + log(w))))
    else
	as.numeric(fu %*% w) # as.*(): drop dimension
}
setMethod("dCopula", signature("matrix",  "mixCopula"), dMixCopula)

## Random Number Generation
setMethod("rCopula", signature("numeric", "mixCopula"),
	  function(n, copula) {
    ## Determine number of samples from each copula
    nj <- as.vector(rmultinom(n = 1, size = n, prob = copula@w))
    ## Draw nj samples from the jth copula
    U <- lapply(seq(along = nj),
                function(j) if(nj[j] > 0) rCopula(nj[j], copula@cops[[j]]) else NULL) # fix to work with do.call() if nj == 0
    ## Bind rows together and randomly permute the entries
    do.call(rbind, U)[sample.int(n), ]
          })

## Tail dependence
setMethod("lambda", "mixCopula", function(copula, ...)
    setNames(c(vapply(copula@cops, lambda, numeric(2)) %*% copula@w),
             c("lower", "upper")))

## Spearman's rho
setMethod("rho", "mixCopula", function(copula, ...) # note: Kendall's tau non-trivial
    c(vapply(copula@cops, rho, 1.1) %*% copula@w))



