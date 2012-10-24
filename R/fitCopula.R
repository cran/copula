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


### Class fitCopula ############################################################

setClass("fitCopula",
	 representation(copula = "copula"),
	 contains="fittedMV" #-> ./Classes.R
	 ## FIXME , validity = function(object) TRUE
	 )

setMethod("paramNames", "fitCopula", function(x) x@copula@param.names)

print.fitCopula <- function(x, digits = max(3, getOption("digits") - 3),
    signif.stars = getOption("show.signif.stars"), ...)
{
    foo <- summary.fitCopula(x)
    cat("fitCopula() estimation based on '", x@method, "'\nand a sample of size ",
	x@nsample, ".\n", sep="")
    printCoefmat(foo$coefficients, digits = digits, signif.stars = signif.stars,
		 na.print = "NA", ...)
    if (!is.na(foo$loglik))
	cat("The maximized loglikelihood is ", format(foo$loglik, digits=digits), "\n")
    if (!is.na(foo$convergence)) {
	if(foo$convergence)
            cat("Convergence problems: code is", foo$convergence, "see ?optim.\n")
	else cat("Optimization converged\n")
    }
    if(!is.null(cnts <- x@fitting.stats$counts) && !all(is.na(cnts))) {
	cat("Number of loglikelihood evaluations:\n"); print(cnts, ...)
    }
    invisible(x)
}

summary.fitCopula <- function(object, ...) {
  estimate <- object@estimate
  se <- sqrt(diag(object@var.est))
  zvalue <- estimate / se
  pvalue <- 2* pnorm(abs(zvalue), lower.tail=FALSE)
  coef <- cbind(Estimate = estimate, "Std. Error" = se,
                "z value" = zvalue, "Pr(>|z|)" = pvalue)
  rownames(coef) <- paramNames(object)
  structure(class = "summary.fitCopula",
            list(method = object@method,
                 loglik = object@loglik,
		 convergence = object@fitting.stats[["convergence"]],
                 coefficients = coef))
}

setMethod("show", signature("fitCopula"),
	  function(object) print.fitCopula(object))


## Wrapper function ############################################################

fitCopula <- function(copula, data, method = c("mpl","ml","itau","irho"),
                      start = NULL, lower = NULL, upper = NULL,
                      optim.method = "BFGS", optim.control = list(maxit=1000),
                      estimate.variance = TRUE, hideWarnings = TRUE)
{
  if(!is.matrix(data)) {
    warning("coercing 'data' to a matrix.")
    data <- as.matrix(data); stopifnot(is.matrix(data))
  }
  switch(match.arg(method),
	 "ml" =
	 fitCopula.ml(copula, data, start=start, lower=lower, upper=upper,
		      method=optim.method, optim.control=optim.control,
		      estimate.variance=estimate.variance, hideWarnings=hideWarnings),
	 "mpl" =
	 fitCopula.mpl(copula, data, start=start, lower=lower, upper=upper,
		       optim.method=optim.method, optim.control=optim.control,
		       estimate.variance=estimate.variance, hideWarnings=hideWarnings),
	 "itau" = fitCopula.itau(copula, data, estimate.variance=estimate.variance),
	 "irho" = fitCopula.irho(copula, data, estimate.variance=estimate.variance))
}


fitCopStart <- function(copula, data, default = copula@parameters)
{
    if (hasMethod("iTau", clc <- class(copula))) {
	ccl <- getClass(clc)
	.par.df <- has.par.df(copula, ccl)
	start <- fitCopula.itau(if(.par.df) as.df.fixed(copula, ccl) else copula,
				data, FALSE, FALSE)@estimate
	if(.par.df) ## add starting value for 'df'
	    start <- c(start, copula@df)
	if(!is.finite(loglikCopula(start, data, copula)))
	    default else start
    }
    else default
}

## fitCopula with maximizing pseudo-likelihood #################################

fitCopula.mpl <- function(copula, u, start=NULL,
                          lower=NULL, upper=NULL,
                          optim.method=NULL, optim.control=list(),
                          estimate.variance = TRUE,
                          hideWarnings = TRUE)
{
    fit <- fitCopula.ml(copula, u, start=start, lower=lower, upper=upper,
                        method=optim.method, optim.control=optim.control,
                        estimate.variance=FALSE, hideWarnings=hideWarnings)
    var.est <- if(estimate.variance)
        varPL(fit@copula, u) / nrow(u) else {
            q <- length(copula@parameters)
            matrix(NA, q, q)
        }
    new("fitCopula",
        estimate = fit@estimate,
        var.est = var.est,
        method = "maximum pseudo-likelihood",
        loglik = fit@loglik,
        fitting.stats = fit@fitting.stats,
        nsample = nrow(u),
        copula = fit@copula)
}


##' fitCopula using inversion of Kendall's tau
fitCopula.itau <- function(copula, x, estimate.variance=TRUE, warn.df=TRUE) {
  ccl <- getClass(class(copula))
  isEll <- extends(ccl, "ellipCopula")
  if(has.par.df(copula, ccl, isEll)) { ## must treat it as "df.fixed=TRUE"
      if(warn.df)
          warning("\"itau\" fitting ==> copula coerced to 'df.fixed=TRUE'")
      copula <- as.df.fixed(copula, classDef = ccl)
  }
  q <- length(copula@parameters)
  tau <- cor(x, method="kendall")
  itau <- fitKendall(copula, tau)
  itau <- itau[lower.tri(itau)]

  X <- getXmat(copula)
  estimate <-
      as.vector(# stripping attributes
		if(isEll && copula@dispstr == "ar1") ## special treatment
		exp(coef(lm(log(itau) ~ X - 1))) else coef(lm(itau ~ X - 1)))
  copula@parameters <- estimate
  var.est <- if (estimate.variance)
    varKendall(copula, x) / nrow(x) else matrix(NA, q, q)
  new("fitCopula",
      estimate = estimate,
      var.est = var.est,
      method = "inversion of Kendall's tau",
      loglik = NA_real_,
      fitting.stats = list(convergence = NA_integer_),
      nsample = nrow(x),
      copula = copula)
}


## fitCopula using inversion of Spearman's rho #################################

fitCopula.irho <- function(copula, x, estimate.variance=TRUE, warn.df=TRUE) {
  ccl <- getClass(class(copula))
  isEll <- extends(ccl, "ellipCopula")
  if(has.par.df(copula, ccl, isEll)) { ## must treat it as "df.fixed=TRUE"
      if(warn.df)
          warning("\"irho\" fitting ==> copula coerced to 'df.fixed=TRUE'")
      copula <- as.df.fixed(copula, classDef = ccl)
  }
  q <- length(copula@parameters)
  rho <- cor(x, method="spearman")
  irho <- fitSpearman(copula, rho)
  irho <- irho[lower.tri(irho)]
  X <- getXmat(copula)
  estimate <-
      as.vector(# stripping attributes
		if (isEll && copula@dispstr == "ar1") ## special treatment
		exp(coef(lm(log(irho) ~ X - 1))) else coef(lm(irho ~ X - 1)))
  copula@parameters <- estimate
  var.est <- if (estimate.variance)
    varSpearman(copula, x)/nrow(x) else matrix(NA, q, q)
  new("fitCopula",
      estimate = estimate,
      var.est = var.est,
      method = "inversion of Spearman's rho",
      loglik = NA_real_,
      fitting.stats = list(convergence = NA_integer_),
      nsample = nrow(x),
      copula = copula)
}


## fitCopula using maximum pseudo-likelihood ###################################

chkParamBounds <- function(copula) {
  d <- length(param <- copula@parameters)
  m <- length(upper <- copula@param.upbnd)
  if(d != m) return(FALSE)
  m <- length(lower <- copula@param.lowbnd)
  if(d != m) return(FALSE)
  ### return TRUE  iff  ("Parameter value out of bound")
  !(any(is.na(param) | param > upper | param < lower))
}

loglikCopula <- function(param, x, copula) {
  stopifnot(length(copula@parameters) == length(param))
  copula@parameters <- param
  if (chkParamBounds(copula))
      sum(dCopula(x, copula, log=TRUE)) else -Inf # was NaN
}

fitCopula.ml <- function(copula, u, start=NULL,
                         lower=NULL, upper=NULL,
                         method=NULL, optim.control=list(),
                         estimate.variance=TRUE,
                         hideWarnings=FALSE,
                         bound.eps = .Machine$double.eps ^ 0.5)
{
  if(any(u < 0) || any(u > 1))
     stop("'u' must be in [0,1] -- probably rather use pobs(.)")
  stopifnot((d <- ncol(u)) >= 2)
  if (copula@dimension != d)
    stop("The dimension of the data and copula do not match")
  if(is.null(start))
    start <- fitCopStart(copula, u)
  if(any(is.na(start))) stop("'start' contains NA values")
  q <- length(copula@parameters)
  if (q != length(start))
    stop(sprintf("The lengths of 'start' (= %d) and copula@parameters (=%d) differ",
		 length(start), q))

  control <- c(optim.control, fnscale = -1)
  control <- control[!vapply(control, is.null, NA)]

  if (!is.null(optim.control[[1]])) control <- c(control, optim.control)
  meth.has.bounds <- method %in% c("Brent","L-BFGS-B")
  if (is.null(lower))
    lower <- if(meth.has.bounds) copula@param.lowbnd + bound.eps else -Inf
  if (is.null(upper))
    upper <- if(meth.has.bounds) copula@param.upbnd  - bound.eps else Inf
  fit <- optim(start, loglikCopula,
               lower=lower, upper=upper,
               method=method,
               copula = copula, x = u,
               control = control)
  copula@parameters[1:q] <- fit$par
  loglik <- fit$val
  if(fit$convergence > 0)
      warning("possible convergence problem: optim() gave code=",
              fit$convergence)
  
  varNA <- matrix(NA_real_, q, q)
  var.est <- if(estimate.variance && fit$convergence == 0) {
    fit.last <- optim(copula@parameters, loglikCopula,
                      lower=lower, upper=upper,
                      method=method,
                      copula=copula, x=u,
                      control=c(control, maxit=0), hessian=TRUE)

    vcov <- tryCatch(solve(-fit.last$hessian), error = function(e) e)
    if(is(vcov, "error")) {
      warning("Hessian matrix not invertible: ", vcov$message)
      varNA
    } else vcov ## ok
  } else varNA

  new("fitCopula",
      estimate = fit$par,
      var.est = var.est,
      method = "maximum likelihood",
      loglik = loglik,
      ## convergence = as.integer(convergence),
      fitting.stats = c(list(method=method),
                     fit[c("convergence", "counts", "message")], control),
      nsample = nrow(u),
      copula = copula)
}


## functions used in estimation and variance computation #######################

## taken from QRMlib and modified
## credit to Alexander McNeil and Scott Ulman
makePosDef <- function (mat, delta = 0.001) {
  decomp <- eigen(mat)
  Lambda <- decomp$values
  if (any(Lambda < 0)) {
    warning("Estimate is not positive-definite. Correction applied.")
    Lambda[Lambda < 0] <- delta
    Gamma <- decomp$vectors
    newmat <- Gamma %*% diag(Lambda) %*% t(Gamma)
    D <- 1/sqrt(diag(newmat))
    diag(D) %*% newmat %*% diag(D)
  }
  else
    mat
}

## rho given as a square matrix
fitSpearman <- function(cop,rho)  {
  stopifnot(is.numeric(p <- ncol(rho)), p == nrow(rho))
  sigma <- matrix(1,p,p)
  for (j in 1:(p-1))
    for (i in (j+1):p)
      {
        sigma[i,j] <- iRho(cop,rho[i,j])
        sigma[j,i] <- sigma[i,j]
      }

  ## make positive definite if necessary
  if (is(cop, "ellipCopula"))
    makePosDef(sigma, delta=0.001)
  else
    sigma
}

## tau given as a square matrix
fitKendall <- function(cop,tau) {
  stopifnot(is.numeric(p <- ncol(tau)), p == nrow(tau))
  sigma <- matrix(1,p,p)
  for (j in 1:(p-1))
    for (i in (j+1):p)
      {
        sigma[i,j] <- iTau(cop,tau[i,j])
        sigma[j,i] <- sigma[i,j]
      }

  ## make positive definite if necessary
  if (is(cop, "ellipCopula"))
    makePosDef(sigma, delta=0.001)
  else
    sigma
}


## variance/covariance of the pseudo-likelihood estimator ######################

influ.terms <- function(u, influ, q)
{
  p <- ncol(u)
  n <- nrow(u)
  ## integral wrt empirical copula (-> sum):
  S <- matrix(0,n,q)
  for (i in 1:p) {
      o.i <- order(u[,i], decreasing=TRUE)
      obi <- rank(u[,i])## FIXME?  add.influ() in  ./gofTests.R  uses  ecdf(.) * M  here
      S <- S + rbind(rep.int(0,q),
                     apply(influ[[i]][o.i,,drop=FALSE],2,cumsum))[n + 1 - obi,,drop=FALSE]
  }
  S / n
}

##' @title variance-covariance (vcov) matrix for Pseudo Likelihood estimat
##' @param cop the *fitted* copula
##' @param u the available pseudo-observations
##' @return vcov matrix
varPL <- function(cop,u)
{
    q <- length(cop@parameters)
    ans <- matrix(NA_real_, q, q)
    ccl <- getClass(clc <- class(cop))
    isEll <- extends(ccl, "ellipCopula")
    ## check if variance can be computed
    msg <- gettext("The variance estimate cannot be computed for this copula.",
                  "  Rather use  'estimate.variance = FALSE'")
    if (!isEll && (!hasMethod("dcopwrap", clc) ||
                   !hasMethod("derPdfWrtArgs", clc))) {
	warning(msg); return(ans)
    }

    p <- cop@dimension
    n <- nrow(u)

    ## influence: second part
    ## integrals computed from the original pseudo-obs u by Monte Carlo
    dcop <- dcopwrap(cop,u) ## wrapper: either dCopula() or '1' (for ellip.)
    ## New    dcop <- dCopula(u,cop) ## in some cases
    ## influ0 <- score (cop,u, dcop)
    ## derArg <- score2(cop,u, dcop)
    influ0 <- derPdfWrtParams(cop,u) / dcop # c. / c  = of dim.  n x q  {== cop<foo>@score()}
    if(is.na(influ0[1])) {
	warning(msg); return(ans)
    }
    derArg <- derPdfWrtArgs(cop,u)   / dcop #         = of dim.  n x p  {missing in copFoo -- TODO}

    ## TODO: use  array instead of list (and change influ.terms() accordingly)
    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    inve <- solve(var(influ0))
    if(has.par.df(cop, ccl, isEll)) {
        ## currently cannot get var/cov for 'df' part
	ans[-q,-q] <- inve %*% var(influ0 - influ.terms(u,influ, q-1)) %*% inve
	ans
    }
    else
	inve %*% var(influ0 - influ.terms(u,influ, q)) %*% inve
}


## variance of the estimator based on Kendall's tau ############################

## copula is the FITTED copula
getL <- function(copula) {
  ## for ellipCopula only
  p <- copula@dimension
  pp <- p * (p - 1) / 2

  dgidx <- outer(1:p, 1:p, "-")
  dgidx <- dgidx[lower.tri(dgidx)]

  if (!is(copula, "ellipCopula")) {
    matrix(1/pp, nrow=pp, ncol=1)
  } else if (copula@dispstr == "ex") {
    matrix(1/pp, nrow=pp, ncol=1)
  } else if(copula@dispstr == "un") {
    diag(pp)
  } else if(copula@dispstr == "toep") {
    mat <- model.matrix(~ factor(dgidx) - 1)
    mat / matrix(colSums(mat), nrow = pp, ncol=p - 1, byrow=TRUE)
  } else if(copula@dispstr == "ar1") {
    stop("Not implemented yet for the dispersion structure 'ar1'.")
    ## estimate log(rho) first and then exponetiate back
    mat <- model.matrix(~ factor(dgidx) - 1)
    mat * matrix(1:(p - 1), nrow=pp, ncol=p - 1, byrow=TRUE)
  }
  else stop("Not implemented yet for the dispersion structure.")
}


getXmat <- function(copula) {
    p <- copula@dimension
    pp <- p * (p - 1) / 2
    if (!is(copula, "ellipCopula")) ## one-parameter non-elliptical copula
	matrix(1, nrow=pp, ncol=1)
    else {
	switch(copula@dispstr,
	       "ex" = matrix(1, nrow=pp, ncol=1),
	       "un" = diag(pp),
	       "toep" =,
	       "ar1" = {
		   dgidx <- outer(1:p, 1:p, "-")
		   dgidx <- dgidx[lower.tri(dgidx)]
		   if(copula@dispstr == "toep")
		       model.matrix(~ factor(dgidx) - 1)
		   else { ## __"ar1"__
		       ## estimate log(rho) first and then exponetiate back
		       ## mat <- model.matrix(~ factor(dgidx) - 1)
		       ## mat %*% diag(1:(p - 1))
		       matrix(dgidx, ncol=1)
		   }
	       },
	       stop("Not implemented yet for this copula/dispersion structure."))
    }
}


varInfluAr1 <- function(cop, v, L, der) {
  ## v is influence for tau or rho
  p <- cop@dimension
  pp <- p * (p - 1) / 2
  n <- nrow(v)

  ## estimate log(r) first, then log(theta), and then exponetiate back
  ## r is the lower.tri of sigma
  sigma <- getSigma(cop) # assuming cop is the fitted copula
  ## influ for log(r)
  r <- sigma[lower.tri(sigma)]
  der <- if (der == "tau") dTauFun(cop)(r) else dRhoFun(cop)(r)
  D <- diag(x = 1 / r / der, pp)
  v <- v %*% D
  ## influ for log(theta)
  v <- v %*% L

  ## influ for theta
  theta <- cop@parameters[1]
  v %*% theta
}

##' Variance of the estimator based on Kendall's tau
##' See Kojadinovic & Yan (2010) Comparison of three semiparametric ... IME 47, 52--63
##' @param cop the \bold{fitted} copula
##' @param u the available pseudo-observations
varKendall <- function(cop,u) {
  ## check if variance can be computed
  if (!hasMethod("dTau", class(cop))) {
    warning("The variance estimate cannot be computed for a copula of class ", class(cop))
    q <- length(cop@parameters)
    return(matrix(NA, q, q))
  }
  p <- cop@dimension
  n <- nrow(u)
  ec <- numeric(n)
  v <- matrix(0,n,p*(p-1)/2)

  ## get influence functions for tau
  l <- 1
  for (j in 1:(p-1)) {
    for (i in (j+1):p) {
      for (k in 1:n) ## can this be vectorized?
        ec[k] <- sum(u[,i] <= u[k,i] & u[,j] <= u[k,j])/n
      v[,l] <- 2 * ec - u[,i] - u[,j]
      l <- l + 1
    }
  }
  ## L <- getL(cop)
  X <- getXmat(cop)
  L <- t(solve(crossprod(X), t(X)))
  v <- if (is(cop, "ellipCopula") && cop@dispstr == "ar1") { ## special treatment
    varInfluAr1(cop, v, L, "tau")
  }
  else {
    ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
    D <- if (length(cop@parameters) == 1) 1 / dTau(cop) else diag(1 / dTau(cop))
    v %*% L %*% D
  }
  16 * var(v)
}


##' Variance of the estimator based on Spearman's rho
##' See Kojadinovic & Yan (2010) Comparison of three semiparametric ... IME 47, 52--63
##'
##' @title variance of the estimator based on Spearman's rho
##' @param cop  the FITTED copula
##' @param u the available pseudo-observations
##' @return
##' @author
varSpearman <- function(cop,u)  {
  ## check if variance can be computed
  if (!hasMethod("dRho", class(cop))) {
    warning("The variance estimate cannot be computed for this copula.")
    q <- length(cop@parameters)
    return(matrix(NA, q, q))
  }

  stopifnot((p <- cop@dimension) >= 2)
  n <- nrow(u)
  v <- matrix(0, n, p*(p-1)/2)

  ord <- apply(u, 2, order, decreasing=TRUE)
  ordb <- apply(u, 2, rank)# ties : "average"
  storage.mode(ordb) <- "integer" # as used below

  l <- 0L
  for (j in 1:(p-1)) {
    for (i in (j+ 1L):p)
      v[,(l <- l + 1L)] <- u[,i] * u[,j] +
          c(0, cumsum(u[ord[,i], j]))[n + 1L - ordb[,i]] / n +
          c(0, cumsum(u[ord[,j], i]))[n + 1L - ordb[,j]] / n
  }
  ## L <- getL(cop)
  X <- getXmat(cop)
  L <- t(solve(crossprod(X), t(X)))
  v <- if (is(cop, "ellipCopula") && cop@dispstr == "ar1") {
      varInfluAr1(cop, v, L, "rho")
  } else {
      ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
      D <- if (length(cop@parameters) == 1) 1 / dRho(cop) else diag(1 / dRho(cop))
      v %*% L %*% D
  }
  144 * var(v)
}

