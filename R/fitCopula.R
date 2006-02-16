loglikCopula <- function(param, x, copula) {
  copula@parameters <- param
  sum(log(dcopula(copula, x)))
}

fitCopula <- function(data, copula, start,
                      optim.control=list(NULL), method="BFGS") {
  if (copula@dimension != ncol(data))
    stop("The dimention of the data and copual not match.\n")
  if (length(copula@parameters) != length(start))
    stop("The length of start and copula parameters not match.\n")

  control <- c(optim.control, fnscale=-1)
  fit <- optim(start, loglikCopula, copula = copula, x = data, control=control)
  if (fit$convergence > 0)
    warning("possible convergence problem: optim gave code=", fit$convergence)
  copula@parameters <- fit$par
  loglik <- fit$val

  fit.last <- optim(copula@parameters, loglikCopula, copula=copula, x =data, control=c(control, maxit=1), hess=TRUE)
    
  ans <- list(est = fit$par,
              var.est = solve(-fit.last$hess),
              loglik = loglik,
              fit = fit)
  ans
}
        
  
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
  sum(log(dmvdc(mvdc, x)))
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
  fit <- optim(start, loglikMvdc, mvdc=mvdc, x = data, control=control)
  if (fit$convergence > 0)
    warning("possible convergence problem: optim gave code=", fit$convergence)
  loglik <- fit$val

  fit.last <- optim(fit$par, loglikMvdc, mvdc=mvdc, x =data, control=c(control, maxit=1), hess=TRUE)
    
  ans <- list(est = fit$par,
              var.est = solve(-fit.last$hess),
              loglik = loglik,
              fit = fit)
  ans
}

