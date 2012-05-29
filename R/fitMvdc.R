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

setClass("fitMvdc",
         representation(estimate = "numeric",
                        var.est = "matrix",
                        loglik = "numeric",
                        convergence = "integer",
                        nsample = "integer",
                        mvdc = "mvdc"),
         validity = function(object) TRUE
         )


setClass("summaryFitMvdc",
         representation(loglik = "numeric",
                        convergence = "integer",
                        parameters = "data.frame"),
         validity = function(object) TRUE
         )

showFitMvdc <- function(object) {
  foo <- summaryFitMvdc(object)
  cat("The Maximum Likelihood estimation is based on ", object@nsample, " observations.\n")
  p <- object@mvdc@copula@dimension
  marNpar <- unlist(lapply(object@mvdc@paramMargins, length))
  idx2 <- cumsum(marNpar)
  idx1 <- idx2 - marNpar + 1
  margid <- object@mvdc@marginsIdentical
  if (sum(marNpar) > 0) { ## sometimes there is no marginal params
    if(margid){
      cat("Margins:\n")
      print(foo@parameters[idx1[1]:idx2[1], 1:4, drop=FALSE])
    }
    else {
      for (i in 1:p) {
        cat("Margin ", i, ":\n")
        print(foo@parameters[idx1[i]:idx2[i], 1:2, drop=FALSE])
      }
    }
  }
  cat("Copula:\n")
  copParIdx <- seq_along(object@mvdc@copula@parameters)
  print(foo@parameters[sum(marNpar) + copParIdx,
                       if(margid) 1:4 else 1:2, drop=FALSE])
  cat("The maximized loglikelihood is ", foo@loglik, "\n")
  cat("The convergence code is ", foo@convergence, "see ?optim.\n")
  invisible(object)
}

summaryFitMvdc <- function(object) {
  estimate <- object@estimate
  se <- sqrt(diag(object@var.est))
  zvalue <- estimate / se
  pvalue <- (1 - pnorm(abs(zvalue))) * 2
  ##ans <- object[c("loglik", "convergence")]
  parameters <- data.frame(estimate, se, zvalue, pvalue)
  marNpar <- unlist(lapply(object@mvdc@paramMargins, length))
  p <- object@mvdc@copula@dimension

  margpnames <-
      if (sum(marNpar) == 0) NULL
      else if(object@mvdc@marginsIdentical)
          c(paste(paste("m", lapply(object@mvdc@paramMargins, names)[[1]], sep=".")))
      else
          c(paste(paste("m", rep(1:p, marNpar), sep=""),
                  unlist(lapply(object@mvdc@paramMargins, names)), sep="."))
  dimnames(parameters) <-
    list(c(margpnames, object@mvdc@copula@param.names),
         c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  new("summaryFitMvdc",
      loglik = object@loglik,
      convergence = object@convergence,
      parameters = parameters)
}


setMethod("show", signature("fitMvdc"), showFitMvdc)
setMethod("summary", signature("fitMvdc"), summaryFitMvdc)


################################################################################

loglikMvdc <- function(param, x, mvdc, suppressMessages=FALSE) {
  p <- mvdc@copula@dimension
  margid <- mvdc@marginsIdentical
  marNpar <- unlist(lapply(mvdc@paramMargins, length))
  idx2 <- cumsum(marNpar)
  idx1 <- idx2 - marNpar + 1

  for (i in 1:p) {
    if (marNpar[i] > 0) {
      ## parnames <- mvdc@paramMargins[[i]]
      k <- if(margid) 1 else i
      par <- param[idx1[k]: idx2[k]]
      ## names(par) <- parnames
      ## mvdc@paramMargins[i] <- as.list(par)
      for (j in 1:marNpar[i]) mvdc@paramMargins[[i]][j] <- par[j]
    }
  }
  mvdc@copula@parameters <-
      if (idx2[p] == 0) # no marginal parameters
          param else if(margid) param[- (1:idx2[1])] else param[- (1:rev(idx2)[1])]

  ## messageOut may be used for debugging
  if (suppressMessages) {
    messageOut <- textConnection("fitMessages", open="w", local=TRUE)
    sink(messageOut); sink(messageOut, type="message")
    options(warn = -1) ## ignore warnings; can be undesirable!
  }

  loglik <- try(sum(log(dmvdc(mvdc, x))))

  if (suppressMessages) {
    options(warn = 0)
    sink(type="message"); sink(); close(messageOut)
  }
  if (inherits(loglik, "try-error")) loglik <- NaN
  loglik
}

fitMvdc <- function(data, mvdc, start,
                    optim.control=list(), method="BFGS")
{
  copula <- mvdc@copula
  if (copula@dimension != ncol(data))
    stop("The dimension of the data and copula do not match.\n")
  marNpar <- unlist(lapply(mvdc@paramMargins, length))
  if(mvdc@marginsIdentical){
    if(length(copula@parameters) + marNpar[1] != length(start))
      stop("The length of start and mvdc parameters do not match.\n")
  }
  else {
    if(length(copula@parameters) + sum(marNpar) != length(start))
      stop("The length of start and mvdc parameters do not match.\n")
  }
  control <- c(optim.control, fnscale=-1)
  notgood <- unlist(lapply(control, is.null))
  control <- control[!notgood]

  fit <- optim(start, loglikMvdc,  method=method, mvdc=mvdc, x = data, suppressMessages=TRUE, control=control)
  if (fit$convergence > 0)
    warning("possible convergence problem: optim gave code=", fit$convergence)
  loglik <- fit$val

  fit.last <- optim(fit$par, loglikMvdc, method=method, mvdc=mvdc, x =data, suppressMessages=TRUE, control=c(control, maxit=1), hessian=TRUE)

  var.est <- try(solve(-fit.last$hessian))
  if (inherits(var.est, "try-error"))
    warning("Hessian matrix not invertible")

  new("fitMvdc",
             estimate = fit$par,
             var.est = var.est,
             loglik = loglik,
             convergence = fit$convergence,
             nsample = nrow(data),
             mvdc = mvdc)
}
