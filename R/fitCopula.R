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
         representation(estimate = "numeric",
                        var.est = "matrix",
                        method = "character",
                        loglik = "numeric",
                        convergence = "integer",
                        optimOpts = "list",
                        nsample = "integer",
                        copula = "copula")
         ## FIXME validity = function(object) TRUE
         )

setClass("summaryFitCopula",
         representation(method = "character",
                        loglik = "numeric",
                        convergence = "integer",
                        parameters = "data.frame")
         ## FIXME validity = function(object) TRUE
         )

showFitCopula <- function(object) {
  foo <- summaryFitCopula(object)
  cat("The estimation is based on the ", object@method, "\nand a sample of size ", object@nsample, ".\n", sep="")
  print(foo@parameters)
  if (!is.na(foo@loglik)) cat("The maximized loglikelihood is ", foo@loglik, "\n")
  if (!is.na(foo@convergence)) cat("The convergence code is ", foo@convergence, "\n")
  invisible(object)
}

summaryFitCopula <- function(object) {
  estimate <- object@estimate
  se <- sqrt(diag(object@var.est))
  zvalue <- estimate / se
  pvalue <- (1 - pnorm(abs(zvalue))) * 2
  parameters <- data.frame(estimate, se, zvalue, pvalue)
  dimnames(parameters) <-
    list(object@copula@param.names,
         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ## MM{FIXME}: S4 class is overkill here
  new("summaryFitCopula",
      method = object@method,
      loglik = object@loglik,
      convergence = object@convergence,
      parameters = parameters)
}

setMethod("show", signature("fitCopula"), showFitCopula)
setMethod("summary", signature("fitCopula"), summaryFitCopula)


## Wrapper function ############################################################

fitCopula <- function(copula, data, method = c("mpl","ml","itau","irho"),
                      start = NULL, lower = NULL, upper = NULL,
                      optim.method = "BFGS", optim.control = list(maxit=1000),
                      estimate.variance = TRUE, hideWarnings = TRUE)
{
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


fitCopStart <- function(copula, data, hideWarnings, default = copula@parameters)
{
  if (hasMethod("calibKendallsTau", clc <- class(copula))) {
    start <- fitCopula.itau(copula, data, FALSE)@estimate
    if(extends(clc, "tCopula") && !copula@df.fixed)
      ## add starting value for 'df'
      start <- c(start, 4)
    if(is.na(loglikCopula(start, data, copula, hideWarnings)))
      default else start
  }
  else default
}

## fitCopula with maximizing pseudo-likelihood #################################

fitCopula.mpl <- function(copula, data, start=NULL,
                          lower=NULL, upper=NULL,
                          optim.method=NULL, optim.control=list(),
                          estimate.variance = TRUE,
                          hideWarnings = TRUE) {
  q <- length(copula@parameters)
  if(is.null(start))
    start <- fitCopStart(copula, data, hideWarnings)

  fit <- fitCopula.ml(copula, data, start=start, lower=lower, upper=upper,
        	      method=optim.method, optim.control=optim.control,
                      estimate.variance=FALSE, hideWarnings=hideWarnings)
  var.est <- if(estimate.variance)
    varPL(fit@copula, data) / nrow(data) else matrix(NA, q, q)
  new("fitCopula",
      estimate = fit@estimate,
      var.est = var.est,
      method = "maximum pseudo-likelihood",
      loglik = fit@loglik,
      convergence = fit@convergence,
      optimOpts = fit@optimOpts,
      nsample = nrow(data),
      copula = fit@copula)
}


## fitCopula using inversion of Kendall's tau ##################################

fitCopula.itau <- function(copula, data, estimate.variance=TRUE) {
  q <- length(copula@parameters)
  X <- getXmat(copula)
  tau <- cor(data, method="kendall")
  itau <- fitKendall(copula, tau)
  itau <- itau[lower.tri(itau)]
  estimate <-
      if (is(copula, "ellipCopula") && copula@dispstr == "ar1") ## special treatment
          exp(coef(lm(log(itau) ~ X - 1))) else coef(lm(itau ~ X - 1))
  attributes(estimate) <- NULL ## strip attributes
  ## Note that '@ parameters' may contain "'df' at end" for tCopula(*, df.fixed=FALSE) :
  copula@parameters[seq_along(estimate)] <- estimate
  var.est <- if (estimate.variance)
    varKendall(copula, data) / nrow(data) else matrix(NA, q, q)
  new("fitCopula",
      estimate = estimate,
      var.est = var.est,
      method = "inversion of Kendall's tau",
      loglik = NA_real_,
      convergence = NA_integer_,
      nsample = nrow(data),
      copula = copula)
}


## fitCopula using inversion of Spearman's rho #################################

fitCopula.irho <- function(copula, data, estimate.variance=TRUE) {
  q <- length(copula@parameters)
  X <- getXmat(copula)
  rho <- cor(data, method="spearman")
  irho <- fitSpearman(copula, rho)
  irho <- irho[lower.tri(irho)]
  estimate <- if (is(copula, "ellipCopula") && copula@dispstr == "ar1") ## special treatment
      exp(coef(lm(log(irho) ~ X - 1))) else coef(lm(irho ~ X - 1))
  attributes(estimate) <- NULL ## strip attributes
  ## Note that '@ parameters' may contain "'df' at end" for tCopula(*, df.fixed=FALSE) :
  copula@parameters[seq_along(estimate)] <- estimate
  var.est <- if (estimate.variance)
    varSpearman(copula, data)/nrow(data) else matrix(NA, q, q)
  new("fitCopula",
      estimate = estimate,
      var.est = var.est,
      method = "inversion of Spearman's rho",
      loglik = as.numeric(NA),
      convergence = as.integer(NA),
      nsample = nrow(data),
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

loglikCopula <- function(param, x, copula, hideWarnings=FALSE) {
  stopifnot(length(copula@parameters) == length(param))
  copula@parameters <- param
  if (!chkParamBounds(copula)) return(NaN)

## FIXME:  use suppressMessages() {and live without  "fitMessages"
  ## messageOut may be used for debugging
  if (hideWarnings) {
    messageOut <- textConnection("fitMessages", open="w", local=TRUE)
    sink(messageOut); sink(messageOut, type="message")
    options(warn = -1) ## ignore warnings; can be undesirable!
  }

  loglik <- try(sum(dcopula(copula, x, log=TRUE)))
  ##old: loglik <- try(sum(log(dcopula(copula, x))))

  if (hideWarnings) {
    options(warn = 0)
    sink(type="message"); sink(); close(messageOut)
  }
  if (inherits(loglik, "try-error")) return(NaN)
  loglik
}

fitCopula.ml <- function(copula, data, start=NULL,
                         lower=NULL, upper=NULL,
                         method=NULL, optim.control=list(),
                         estimate.variance=TRUE,
                         hideWarnings=FALSE)
{
  if (copula@dimension != ncol(data))
    stop("The dimension of the data and copula do not match.\n")
  if(is.null(start))
    start <- fitCopStart(copula, data, hideWarnings)
  q <- length(copula@parameters)
  if (q != length(start))
    stop(sprintf("The lengths of 'start' (= %d) and copula@parameters (=%d) differ",
		 length(start), q))

  control <- c(optim.control, fnscale = -1)
  control <- control[!vapply(control, is.null, NA)]

  if (!is.null(optim.control[[1]])) control <- c(control, optim.control)
  eps <- .Machine$double.eps ^ 0.5
  meth.has.bounds <- method %in% c("Brent","L-BFGS-B")
  if (is.null(lower))
    lower <- if(meth.has.bounds) copula@param.lowbnd + eps else -Inf
  if (is.null(upper))
    upper <- if(meth.has.bounds) copula@param.upbnd  - eps else Inf
##  if (p >= 2) {
  if (TRUE) {
    fit <- optim(start, loglikCopula,
                 lower=lower, upper=upper,
                 method=method,
                 copula = copula, x = data, hideWarnings = hideWarnings,
                 control = control)
    copula@parameters[1:q] <- fit$par
    loglik <- fit$val
    convergence <- fit$convergence
    if(convergence > 0)
      warning("possible convergence problem: optim gave code=", convergence)
  }
##   else {
##     fit <- optimize(loglikCopula, lower=lower, upper=upper, maximum=TRUE,
##                     x=data, copula=copula)
##     copula@parameters <- fit$maximum
##     loglik <- fit$objective
##     convergence <- 0
##   }

  if (estimate.variance){
    fit.last <- optim(copula@parameters, loglikCopula,
                      lower=lower, upper=upper,
                      method=method,
                      copula=copula, x=data, hideWarnings=hideWarnings,
                      control=c(control, maxit=0), hessian=TRUE)
    var.est <- try(solve(-fit.last$hessian))
    if (inherits(var.est, "try-error"))
      warning("Hessian matrix not invertible")
  } else var.est <- matrix(NA, q, q)

  new("fitCopula",
      estimate = fit$par,
      var.est = var.est,
      method = "maximum likelihood",
      loglik = loglik,
      convergence = as.integer(convergence),
      optimOpts = list(method = method, control = control),
      nsample = nrow(data),
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
        sigma[i,j] <- calibSpearmansRho(cop,rho[i,j])
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
        sigma[i,j] <- calibKendallsTau(cop,tau[i,j])
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

  o <- ob <- matrix(0,n,p)
  for (i in 1:p)
    {
      o[,i] <- order(u[,i], decreasing=TRUE)
      ob[,i] <- rank(u[,i])
    }

  out <- matrix(0,n,q)
  for (i in 1:p)
      out <- out + rbind(rep(0,q),apply(influ[[i]][o[,i],,drop=FALSE],2,cumsum))[n + 1 - ob[,i],,drop=FALSE] / n
  return(out)
}

## cop is the FITTED copula
## u are the available pseudo-observations
varPL <- function(cop,u)
  {
    ## check if variance can be computed
    if (!hasMethod("dcopwrap", class(cop))) {
      warning("The variance estimate cannot be computed for this copula.")
      q <- length(cop@parameters)
      return(matrix(NA, q, q))
    }

    p <- cop@dimension
    n <- nrow(u)

    ## influence: second part
    ## integrals computed from the original pseudo-obs u by Monte Carlo
    dcop <- dcopwrap(cop,u) ## wrapper
    influ0 <- derPdfWrtParams(cop,u)/dcop
    derArg <- derPdfWrtArgs(cop,u)/dcop

    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    q <- length(cop@parameters)
    inve <- solve(var(influ0))

    return(inve %*% var(influ0 - influ.terms(u,influ,q)) %*% inve)
  }


## variance of the estimator based on Kendall's tau ############################

## cop is the FITTED copula
## u are the available pseudo-observations
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

  dgidx <- outer(1:p, 1:p, "-")
  dgidx <- dgidx[lower.tri(dgidx)]

  if (!is(copula, "ellipCopula")) { ## one-parameter non-elliptical copula
    matrix(1, nrow=pp, ncol=1)
  } else if (copula@dispstr == "ex") {
    matrix(1, nrow=pp, ncol=1)
  } else if(copula@dispstr == "un") {
    diag(pp)
  } else if(copula@dispstr == "toep") {
    model.matrix(~ factor(dgidx) - 1)
  } else if(copula@dispstr == "ar1") {
    ## estimate log(rho) first and then exponetiate back
    ## mat <- model.matrix(~ factor(dgidx) - 1)
    ## mat %*% diag(1:(p - 1))
    matrix(dgidx, ncol=1)
  }
  else stop("Not implemented yet for this copula/dispersion structure.")
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
  der <- if (der == "tau") tauDerFun(cop)(r) else rhoDerFun(cop)(r)
  D <- diag(x = 1 / r / der, pp)
  v <- v %*% D
  ## influ for log(theta)
  v <- v %*% L

  ## influ for theta
  theta <- cop@parameters[1]
  v %*% theta
}

varKendall <- function(cop,u) {
  ## check if variance can be computed
  if (!hasMethod("tauDer", class(cop))) {
    warning("The variance estimate cannot be computed for this copula.")
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
  ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
  D <- if (length(cop@parameters) == 1) 1 / tauDer(cop) else diag(1 / tauDer(cop))
  if (is(cop, "ellipCopula") && cop@dispstr == "ar1") { ## special treatment
    v <- varInfluAr1(cop, v, L, "tau")
    return(16 * var(v))
  } else return(16 * var(v %*% L %*% D))
}


## variance of the estimator based on Spearman's rho ###########################

## cop is the FITTED copula
## u are the available pseudo-observations
varSpearman <- function(cop,u)  {
  ## check if variance can be computed
  if (!hasMethod("rhoDer", class(cop))) {
    warning("The variance estimate cannot be computed for this copula.")
    q <- length(cop@parameters)
    return(matrix(NA, q, q))
  }

  stopifnot((p <- cop@dimension) >= 2)
  n <- nrow(u)
  v <- matrix(0,n,p*(p-1)/2)

  ord <- apply(u, 2, order, decreasing=TRUE)
  ordb <- apply(u, 2, rank)

  l <- 1
  for (j in 1:(p-1)) {
    for (i in (j+1):p)  {
      v[,l] <- u[,i] * u[,j] +
          c(0, cumsum(u[ord[,i], j]))[n + 1 - ordb[,i]] / n +
          c(0, cumsum(u[ord[,j], i]))[n + 1 - ordb[,j]] / n
      l <- l + 1
    }
  }
  ## L <- getL(cop)
  X <- getXmat(cop)
  L <- t(solve(crossprod(X), t(X)))
  if (is(cop, "ellipCopula") && cop@dispstr == "ar1") {
      v <- varInfluAr1(cop, v, L, "rho")
      144 * var(v)
  } else {
      ## Caution: diag(0.5) is not a 1x1 matrix of 0.5!!! check it out.
      D <- if (length(cop@parameters) == 1) 1 / rhoDer(cop) else diag(1 / rhoDer(cop))
      144 * var(v %*% L %*% D)
  }
}

