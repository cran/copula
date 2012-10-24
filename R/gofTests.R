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

##' Goodness-of-fit test for copulas: wrapper function
##' Calls either gofMCLT.PL, gofMCLT.KS or gofPB
##'
##' @title Goodness-of-fit test for copulas: wrapper function
##' @param copula is a copula of the desired family
##' @param x the data
##' @param N the number of bootstrap or multiplier replications
##' @param estim.method estimation method for the unknown parameter
##' @param simulation parametric bootstrap or multiplier
##' @param verbose display progress bar if TRUE
##' @param print.every is deprecated
##' @param optim.method for fitting
##' @param optim.control for fitting
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
gofCopula <- function(copula, x, N = 1000,
                      method = c("Sn", "AnChisq", "AnGamma", "SnB", "SnC"),
                      estim.method = c("mpl", "ml", "irho", "itau"),
                      simulation = c("pb", "mult"),
		      verbose = TRUE, print.every = NULL,
                      optim.method = "BFGS", optim.control = list(maxit=20))
{
    if(!is.matrix(x)) {
        warning("coercing 'x' to a matrix")
        x <- as.matrix(x); stopifnot(is.matrix(x))
    }
    n <- nrow(x)
    p <- ncol(x)

    if (p < 2) stop("The data should be at least of dimension 2")
    if (n < 2) stop("There should be at least 2 observations")

    if (copula@dimension != p)
        stop("The copula and the data should be of the same dimension")

    if (!is.null(print.every)) {
        warning("Argument 'print.every' is deprecated. Please use 'verbose' instead.")
        verbose <- print.every > 0
    }
    method <- match.arg(method)
    estim.method <- match.arg(estim.method)

    gof <-
        switch(match.arg(simulation),
               "pb" = { ## parametric bootstrap
		   gofPB(copula, x, N=N, method=method, estim.method=estim.method,
			 verbose=verbose,
			 optim.method=optim.method, optim.control=optim.control)
               },
               "mult" = { ## multiplier
                   if(method != "Sn")
                       stop("Multiplier method for ", method, " not yet implemented")
                   if (estim.method == "mpl")
                       gofMCLT.PL(copula, x, N=N,
                                  optim.method=optim.method, optim.control=optim.control)
                   else if (estim.method %in% c("irho","itau")) {
                       if (copula@dimension != 2)
                           stop("The simulation method 'mult' can be used in combination with the estimation methods 'irho' and 'itau' only in the bivariate case.")
                       gofMCLT.KS(copula, x, N=N, estim.method=estim.method)
                   }
                   else
                       stop(sprintf("Invalid estimation method '%s'", estim.method))
               },
               ## otherwise:
               stop("Invalid simulation method ", match.arg(simulation)))
    gof
}

##' @title Compute Sn (Cramer-von Mises) statistic; = (2) in Genest et al. (2009)
##' @param u (n x d) matrix of copula observations
##' @param copula 
##' @return Sn
##' @author Martin Maechler
SnCvM <- function(u, copula) {
    .C(cramer_vonMises,
       n = as.integer(nrow(u)),
       p = as.integer(ncol(u)),
       as.double(u),
       as.double(pCopula(u, copula)),
       stat = double(1))$stat
}

##' Goodness-of-fit test based on the parametric bootstrap
##' as proposed by Genest et al. (2008)
##'
##' @title Goodness-of-fit test based on the parametric bootstrap
##' @param copula is a copula of the desired family whose parameters,
##' if necessary, will be used as starting values in fitCopula
##' @param x the data
##' @param N the number of bootstrap replications
##' @param estim.method estimation method for the unknown parameter
##' @param verbose display progress bar is TRUE
##' @param optim.method for fitting
##' @param optim.control for fitting
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
gofPB <- function(copula, x, N, method, estim.method, verbose, optim.method, optim.control)
{
    n <- nrow(x)
    p <- ncol(x)
    ## make pseudo-observations
    u <- pobs(x)

    ## fit the copula and compute the test statistic
    fcop <- fitCopula(copula, u, method=estim.method, estimate.variance=FALSE,
                      optim.method=optim.method, optim.control=optim.control)@copula
    
    T <- if(method == "Sn") {
        SnCvM(u, fcop)
    } else {
        ## transform + test statistic
        gofTstat(rtrafo(u, fcop), method=method)
    }

    ## simulation of the null distribution
    T0 <- numeric(N)
    if (verbose) {
	pb <- txtProgressBar(max = N, style = if(isatty(stdout())) 3 else 1) # setup progress bar
	on.exit(close(pb)) # and close it on exit
    }
    # if (print.every > 0)
    #    cat(paste("Progress will be displayed every", print.every, "iterations.\n"))
    for (i in 1:N) {
        ## (4.1)
        u0 <- pobs(rCopula(n, fcop))

        ## fit copula and compute S_n
	fcop0 <- fitCopula(copula, u0, method=estim.method, estimate.variance=FALSE,
			   optim.method=optim.method, optim.control=optim.control)@copula

        T0[i] <- if(method == "Sn") {
            SnCvM(u0, fcop0)
        } else {
            ## transform + test statistic
            gofTstat(rtrafo(u0, fcop0), method=method)
        }

        if (verbose) setTxtProgressBar(pb, i) # update progress bar
    }

    structure(class = "htest",
	      list(method = sprintf(
		   "Parametric bootstrap based GOF test with 'method'=\"%s\", 'estim.method'=\"%s\"",
		   method, estim.method),
                   parameter = c(parameter = fcop@parameters),
                   statistic = c(statistic = T),
                   p.value = (sum(T0 >= T) +0.5)/(N +1),
                   data.name = deparse(substitute(x))))
}


### Utility function: additional influence terms
influ.add <- function(x0, y0, influ1, influ2)
{
  M <- nrow(y0)
  o1 <- order(y0[,1], decreasing=TRUE)
  o1b <- ecdf(y0[,1])(x0[,1]) * M
  o2 <- order(y0[,2], decreasing=TRUE)
  o2b <- ecdf(y0[,2])(x0[,2]) * M
  return(c(0,cumsum(influ1[o1]))[M + 1 - o1b] / M - mean(influ1 * y0[,1]) +
         c(0,cumsum(influ2[o2]))[M + 1 - o2b] / M - mean(influ2 * y0[,2]))
}

##' Goodness-of-fit test based on the multiplier approach
##' and rank correlation coefficients
##'
##' @title Multiplier GOF with rank correlation coefficients
##' @param cop is a copula of the desired family
##' @param x the data
##' @param N number of multiplier replications
##' @param estim.method fitting method
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
gofMCLT.KS <- function(cop, x, N, estim.method)
  {
    stopifnot(estim.method %in% c("irho","itau"), is.matrix(x))
    n <- nrow(x)
    p <- 2

    ## make pseudo-observations
    u <- pobs(x)

    ## fit the copula {*not* calling optim() ..}
    cop <- fitCopula(cop, u, method=estim.method, estimate.variance=FALSE)@copula

    ## compute the test statistic
    s <- .C(cramer_vonMises_grid,
            as.integer(p),
            as.double(u),
            as.integer(n),
            as.double(u),
            as.integer(n),
            as.double(pCopula(u, cop)),
            stat = double(1))$stat

    ## prepare influence coefficients
    if (estim.method == "itau") ## kendall's tau
        influ <- 4 * (2 * pCopula(u, cop) - u[,1] - u[,2]
                      + (1 - tau(cop))/2) / dTau(cop)
    else if (estim.method == "irho") ## Spearman's rho
        influ <- (12 * (u[,1] * u[,2] + influ.add(u, u, u[,2],u[,1])) -
                  3 - rho(cop)) / dRho(cop)

    ## simulate under H0
    s0 <- .C(multiplier,
             as.integer(p),
             as.double(u),
             as.integer(n),
             as.double(u),
             as.integer(n),
             as.double(dCdtheta(cop,u) %*% influ),
             as.integer(N),
             s0 = double(N))$s0

    structure(class = "htest",
              list(method = sprintf("Multiplier GOF test with 'estim.method'=\"%s\"",
                   estim.method),
                   parameter = c(parameter = cop@parameters),
                   statistic = c(statistic = s),
                   p.value = (sum(s0 >= s)+0.5)/(N+1),
                   data.name = deparse(substitute(x))))
}

### Utility function: Influence coefficients for the multipler gof based on MPL
influCoef <- function(cop,u,v)
{
    p <- cop@dimension

    ## influence: second part
    ## integrals computed from M realizations by Monte Carlo
    M <- nrow(v)
    dcop <- dcopwrap(cop,v) ## wrapper
    influ0 <- derPdfWrtParams(cop,v)/dcop
    derArg <- derPdfWrtArgs  (cop,v)/dcop

    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    ## expectation  e := crossprod(influ0)/M
    solve(crossprod(influ0)/M, t(derPdfWrtParams(cop,u)/dcopwrap(cop,u) -
                                 add.influ(u,v, influ=influ, q = length(cop@parameters))))
}

### Utility function: Second part of influence coefficients
add.influ <- function(u, v, influ, q)
{
  M <- nrow(v)
  p <- ncol(v)
  n <- nrow(u)
  S <- matrix(0,n,q)
  for (i in 1:p) {
      vi <- v[,i]
      o.i <- order(vi, decreasing=TRUE)
      obi <- ecdf(vi)(u[,i]) * M # "FIXME": use findInterval(); keep obi integer throughout
      S <- S + rbind(rep.int(0,q),
                     apply(influ[[i]][o.i,,drop=FALSE],2,cumsum))[M + 1 - obi,,drop=FALSE] / M -
                         matrix(colMeans(influ[[i]] * vi), n,q, byrow=TRUE)
	#matrix(apply(influ[[i]] * vi,2,mean),n,q,byrow=TRUE)
  }
  S
}

##' Multiplier GOF based on MPL
##'
##' @title Multiplier GOF based on MPL
##' @param cop is a copula of the desired family
##' @param x the data
##' @param N number of multiplier replications
##' @param optim.method for MPL fitting
##' @param optim.control for MPL fitting
##' @return an object of class 'htest'
##' @author Ivan Kojadinovic
gofMCLT.PL <- function(cop, x, N, optim.method, optim.control)
{
    n <- nrow(x)
    p <- ncol(x)

    ## make pseudo-observations
    u <- pobs(x)

    ## fit the copula
    cop <- fitCopula(cop, u, method="mpl", estimate.variance=FALSE,
                     optim.method=optim.method, optim.control=optim.control)@copula

    ## compute the test statistic
    s <- .C(cramer_vonMises_grid,
            as.integer(p),
            as.double(u),
            as.integer(n),
            as.double(u),
            as.integer(n),
            as.double(pCopula(u, cop)),
            stat = double(1))$stat

    ## simulate under H0
    s0 <- .C(multiplier,
             as.integer(p),
             as.double(u),
             as.integer(n),
             as.double(u),
             as.integer(n),
             as.double(dCdtheta(cop,u) %*% influCoef(cop,u,u)),
             as.integer(N),
             s0 = double(N))$s0

    structure(class = "htest",
              list(method = "Multiplier GOF test with 'method'=\"mpl\"",
                   parameter = c(parameter = cop@parameters),
                   statistic = c(statistic = s),
                   p.value = (sum(s0 >= s)+0.5)/(N+1),
                   data.name = deparse(substitute(x))))
}


