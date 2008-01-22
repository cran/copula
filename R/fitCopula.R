#################################################################################
##
##   R package Copula by Jun Yan Copyright (C) 2008
##
##   This file is part of the R package copula.
##
##   The R package copula is free software: you can redistribute it and/or modify
##   it under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package copula is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package copula. If not, see <http://www.gnu.org/licenses/>.
##
#################################################################################


setClass("fitCopula",
         representation(est = "numeric",
                        var.est = "matrix",
                        loglik = "numeric",
                        convergence = "integer",
                        nsample = "integer",
                        copula = "copula"),
         validity = function(object) TRUE,
         contains = list()
         )

loglikCopula <- function(param, x, copula) {
  copula@parameters <- param
  loglik <- try(sum(log(dcopula(copula, x))))
  if (inherits(loglik, "try-error")) loglik <- NaN
  loglik
}

fitCopula <- function(data, copula, start,
                      lower=NULL, upper=NULL,
                      optim.control=list(NULL),
                      method="BFGS") {
  if (copula@dimension != ncol(data))
    stop("The dimention of the data and copual not match.\n")
  if (length(copula@parameters) != length(start))
    stop("The length of start and copula parameters not match.\n")

  control <- c(optim.control, fnscale=-1)
  p <- length(copula@parameters)
  eps <- .Machine$double.eps ^ 0.5
  if (is.null(lower)) lower <- ifelse(method == "L-BFGS-B", copula@param.lowbnd + eps, -Inf)
  if (is.null(upper)) upper <- ifelse(method == "L-BFGS-B", copula@param.upbnd - eps, Inf)
##  if (p >= 2) {
  if (TRUE) {
    fit <- optim(start, loglikCopula,
                 lower=lower, upper=upper,
                 method=method,
                 copula = copula, x = data,
                 control=control)
    copula@parameters <- fit$par
    loglik <- fit$val
    convergence <- fit$convergence
    if (fit$convergence > 0)
      warning("possible convergence problem: optim gave code=", fit$convergence)
  }
##   else {
##     fit <- optimize(loglikCopula, lower=lower, upper=upper, maximum=TRUE,
##                     x=data, copula=copula)
##     copula@parameters <- fit$maximum
##     loglik <- fit$objective
##     convergence <- 0
##   }


  fit.last <- optim(copula@parameters, loglikCopula,
                    lower=lower, upper=upper,
                    method=method,
                    copula=copula, x =data,
                    control=c(control, maxit=0), hessian=TRUE)
  var.est <- try(solve(-fit.last$hessian))
  if (inherits(var.est, "try-error"))
    warning("Hessian matrix not invertible")
  ans <- new("fitCopula",
             est = fit.last$par,
             var.est = var.est,
             loglik = loglik,
             convergence = as.integer(convergence),
             nsample = nrow(data),
             copula = copula)
  ans
}
      
  

setClass("fitMvdc",
         representation(est = "numeric",
                        var.est = "matrix",
                        loglik = "numeric",
                        convergence = "integer",
                        nsample = "integer",
                        mvdc = "mvdc"),
         validity = function(object) TRUE,
         contains = list()         
         )



loglikMvdc <- function(param, x, mvdc) {
  p <- mvdc@copula@dimension
  marNpar <- unlist(lapply(mvdc@paramMargins, length))
  idx2 <- cumsum(marNpar)
  idx1 <- idx2 - marNpar + 1
  for (i in 1:p) {
    if (marNpar[i] > 0) {
      ## parnames <- mvdc@paramMargins[[i]]
      par <- param[idx1[i]: idx2[i]]
      ## names(par) <- parnames
      ## mvdc@paramMargins[i] <- as.list(par)
      for (j in 1:marNpar[i]) mvdc@paramMargins[[i]][j] <- par[j]
    }      
  }
  mvdc@copula@parameters <- param[- (1:rev(idx2)[1])]
  loglik <- try(sum(log(dmvdc(mvdc, x))))
  if (inherits(loglik, "try-error")) loglik <- NaN
  loglik
}

fitMvdc <- function(data, mvdc, start,
                    optim.control=list(NULL), method="BFGS") {
  copula <- mvdc@copula
  if (copula@dimension != ncol(data))
    stop("The dimention of the data and copual not match.\n")
  marNpar <- unlist(lapply(mvdc@paramMargins, length))
  if (length(copula@parameters) + sum(marNpar) != length(start))
    stop("The length of start and mvdc parameters not match.\n")

  control <- c(optim.control, fnscale=-1)
  fit <- optim(start, loglikMvdc,  method=method, mvdc=mvdc, x = data, control=control)
  if (fit$convergence > 0)
    warning("possible convergence problem: optim gave code=", fit$convergence)
  loglik <- fit$val

  fit.last <- optim(fit$par, loglikMvdc, method=method, mvdc=mvdc, x =data, control=c(control, maxit=1), hessian=TRUE)
    
  var.est <- try(solve(-fit.last$hessian))
  if (inherits(var.est, "try-error"))
    warning("Hessian matrix not invertible")

  ans <- new("fitMvdc",
             est = fit$par,
             var.est = var.est,
             loglik = loglik,
             convergence = fit$convergence,
             nsample = nrow(data),
             mvdc = mvdc)
  ans
}

setClass("summaryFitCopula",
         representation(loglik = "numeric",
                        convergence = "integer",
                        parameters = "data.frame"),
         validity = function(object) TRUE,
         contains = list()
         )

summaryFitCopula <- function(object) {
  est <- object@est
  se <- sqrt(diag(object@var.est))
  zvalue <- est / se
  pvalue <- (1 - pnorm(abs(zvalue))) * 2
  parameters <- data.frame(est, se, zvalue, pvalue)
  dimnames(parameters) <-
    list(object@copula@param.names,
         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ret <- new("summaryFitCopula",
             loglik = object@loglik,
             convergence = object@convergence,
             parameters = parameters)
  ret
}

showFitCopula <- function(object) {
  foo <- summaryFitCopula(object)
  cat("The ML estimation is based on ", object@nsample, " observations.\n")
  print(foo@parameters)
  cat("The maximized loglikelihood is ", foo@loglik, "\n")
  cat("The convergence code is ", foo@convergence, "\n")
}


setClass("summaryFitMvdc",
         representation(loglik = "numeric",
                        convergence = "integer",
                        parameters = "data.frame"),
         validity = function(object) TRUE,
         contains = list()
         )

summaryFitMvdc <- function(object) {
  est <- object@est
  se <- sqrt(diag(object@var.est))
  zvalue <- est / se
  pvalue <- (1 - pnorm(abs(zvalue))) * 2
  ##ans <- object[c("loglik", "convergence")]
  parameters <- data.frame(est, se, zvalue, pvalue)
  marNpar <- unlist(lapply(object@mvdc@paramMargins, length))
  p <- object@mvdc@copula@dimension

  pnames <- c(paste(paste("m", rep(1:p, marNpar), sep=""),
                    unlist(lapply(object@mvdc@paramMargins, names)), sep="."),
              object@mvdc@copula@param.names)
  
  dimnames(parameters) <-
    list(pnames,
         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  ret <- new("summaryFitMvdc",
             loglik = object@loglik,
             convergence = object@convergence,
             parameters = parameters)
  ret
}

showFitMvdc <- function(object) {
  foo <- summaryFitMvdc(object)
  cat("The ML estimation is based on ", object@nsample, " observations.\n")
  p <- object@mvdc@copula@dimension
  marNpar <- unlist(lapply(object@mvdc@paramMargins, length))
  idx2 <- cumsum(marNpar)
  idx1 <- idx2 - marNpar + 1
  for (i in 1:p) {
    cat("Margin ", i, ":\n")
    print(foo@parameters[idx1[i]:idx2[i], 1:2, drop=FALSE])
  }
  cat("Copula:\n")
  print(foo@parameters[- (1:rev(idx2)[1]), 1:2, drop=FALSE])
  
  cat("The maximized loglikelihood is ", foo@loglik, "\n")
  cat("The convergence code is ", foo@convergence, "\n")
}

setMethod("show", signature("fitCopula"), showFitCopula)
setMethod("show", signature("fitMvdc"), showFitMvdc)

setMethod("summary", signature("fitCopula"), summaryFitCopula)
setMethod("summary", signature("fitMvdc"), summaryFitMvdc)
          
