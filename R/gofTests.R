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


### Goodness-of-fit test for copulas: wrapper function
### Calls either gofMCLT.PL, gofMCLT.KS or gofPB

## copula is a copula of the desired family

gofCopula <- function(copula, x, N = 1000, method = "mpl",
                      simulation = c("pb", "mult"), print.every = 100,
                      optim.method = "BFGS", optim.control = list(maxit=20))
{
    M <- -1 ## fixed - for gofMCLT.*()

    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)

    if (p < 2) stop("The data should be at least of dimension 2")
    if (n < 2) stop("There should be at least 2 observations")

    if (copula@dimension != p)
      stop("The copula and the data should be of the same dimension")

    gof <-
        switch(match.arg(simulation),
               "pb" = { ## parametric bootstrap
                   gofPB(copula, x, N=N, method = method,
                         print.every=print.every, optim.method=optim.method, optim.control=optim.control)
               },
               "mult" = { ## multiplier
                   if (method == "mpl")
                       gofMCLT.PL(copula, x, N=N, M=M,
                                  optim.method=optim.method, optim.control=optim.control)
                   else if (method %in% c("irho","itau")) {
                       if (copula@dimension != 2)
                           stop("The simulation method 'mult' can be used in combination with the estimation methods 'irho' and 'itau' only in the bivariate case.")
                       gofMCLT.KS(copula, x, N=N, method=method, M=M)
                   }
                   else
                       stop(sprintf("Invalid estimation method '%s'", method))
               },
               ## otherwise:
               stop("Invalid simulation method ", match.arg(simulation)))
    class(gof) <- "gofCopula"
    gof
}

print.gofCopula <- function(x, ...)
{
  cat("\nParameter estimate(s):", x$parameters, "\n")
  cat("Cramer-von Mises statistic:", x$statistic,
      "with p-value", x$pvalue, "\n\n")
}


### Goodness-of-fit test based on the parametric bootstrap
### as proposed by Genest et al. (2008)

## copula is a copula of the desired family whose parameters, if necessary,
## will be used as starting values in fitCopula

gofPB <- function(copula, x, N, method, print.every, optim.method, optim.control)
  {
    n <- nrow(x)
    p <- ncol(x)

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    fcop <- fitCopula(copula, u, method, estimate.variance=FALSE,
                      optim.method=optim.method, optim.control=optim.control)@copula

    ## compute the test statistic
    s <- .C(cramer_vonMises,
            as.integer(n),
            as.integer(p),
            as.double(u),
            as.double(pcopula(fcop,u)),
            stat = double(1))$stat

    ## simulation of the null distribution
    s0 <- numeric(N)
    if (print.every > 0)
      cat(paste("Progress will be displayed every", print.every, "iterations.\n"))
    for (i in 1:N)
      {
        if (print.every > 0 && i %% print.every == 0)
          cat(paste("Iteration",i,"\n"))
        u0 <- apply(rcopula(fcop,n),2,rank)/(n+1)

        ## fit the copula
        fcop0 <-  fitCopula(copula, u0, method, estimate.variance=FALSE,
                            optim.method=optim.method, optim.control=optim.control)@copula

        s0[i] <- .C(cramer_vonMises,
                    as.integer(n),
                    as.integer(p),
                    as.double(u0),
                    as.double(pcopula(fcop0,u0)),
                    stat = double(1))$stat
      }

    return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1),
                parameters=fcop@parameters))
  }


### Goodness-of-fit test based on the multiplier approach
### and rank correlation coefficients

## additional influence terms ##################################################

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

gofMCLT.KS <- function(cop, x, N, method, M)
  {
    stopifnot(method %in% c("irho","itau"), is.matrix(x))
    n <- m <- nrow(x)
    p <- 2

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula {*not* calling optim() ..}
    cop <- fitCopula(cop, u, method=method, estimate.variance=FALSE)@copula

    ## grid points where to evaluate the process
    g <- u # pseudo-observations
    G <- n

    pcop <- pcopula(cop, g)

    ## compute the test statistic
    stat <- .C(cramer_vonMises_2,
               as.integer(p),
               as.double(u),
               as.integer(n),
               as.double(g),
               as.integer(G),
               as.double(pcop),
               stat = double(1))$stat

    x0 <-  u # rcopula(cop,m)

    ## prepare influence coefficients
    if (method == "itau") ## kendall's tau
      influ <- 4 * (2 * pcopula(cop,x0) - x0[,1] - x0[,2] + (1 - kendallsTau(cop))/2) / tauDer(cop)
    else if (method == "irho") ## Spearman's rho
      {
        ## integrals computed from M realizations by Monte Carlo
        y0 <- if (M > 0) rcopula(cop,M) else u
        influ <- (12 * (x0[,1] * x0[,2] + influ.add(x0, y0, y0[,2],y0[,1])) -
                  3 - spearmansRho(cop)) / rhoDer(cop)
      }

    ## Simulate under H0
    s0 <- .C(multiplier,
             as.integer(p),
             as.double(x0),
             as.integer(m),
             as.double(g),
             as.integer(G),
             as.double(derCdfWrtParams(cop,g) %*% influ),
             as.integer(N),
             s0 = double(N))$s0

    list(statistic=stat, pvalue=(sum(s0 >= stat)+0.5)/(N+1),
                parameters=cop@parameters)
}


## Multivariate multipler gof based on MPL #####################################

## influence coefficients

influCoef <- function(cop,u,v)
{
    p <- cop@dimension

    ## influence: second part
    ## integrals computed from M realizations by Monte Carlo
    M <- nrow(v)
    dcop <- dcopwrap(cop,v) ## wrapper
    influ0 <- derPdfWrtParams(cop,v)/dcop
    derArg <- derPdfWrtArgs(cop,v)/dcop

    influ <- vector("list",p)
    for (i in 1:p)
        influ[[i]] <- influ0 * derArg[,i]

    ## expectation
    q <- length(cop@parameters)
    e <- crossprod(influ0)
    e <- e/M

    ## MM: FIXME --- *after* we have testing code in ./tests/ !!!
    return(solve(e) %*% t(derPdfWrtParams(cop,u)/dcopwrap(cop,u) - add.influ(u,v,influ,q)))
}


## second part of influence coefficients

add.influ <- function(u, v, influ, q)
{
  M <- nrow(v)
  p <- ncol(v)
  n <- nrow(u)

  o <- matrix(0,M,p)
  ob <- matrix(0,n,p)
  for (i in 1:p)
    {
      o[,i] <- order(v[,i], decreasing=TRUE)
      ob[,i] <- ecdf(v[,i])(u[,i]) * M
    }

  out <- matrix(0,n,q)
  for (i in 1:p)
      out <- out + rbind(rep(0,q),apply(influ[[i]][o[,i],,drop=FALSE],2,cumsum))[M + 1 - ob[,i],,drop=FALSE] / M -
        matrix(colMeans(influ[[i]] * v[,i]),n,q,byrow=TRUE)
        #matrix(apply(influ[[i]] * v[,i],2,mean),n,q,byrow=TRUE)
  return(out)
}

## goodness-of-fit test

## cop is a copula of the desired family whose parameters, if necessary, will be used
## as starting values in fitCopula

gofMCLT.PL <- function(cop, x, N, M, optim.method, optim.control)
{
    n <- m <- nrow(x)
    p <- ncol(x)

    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    cop <- fitCopula(cop, u, method="mpl", estimate.variance=FALSE,
                     optim.method=optim.method, optim.control=optim.control)@copula

    ## grid points where to evaluate the process
    g <- u  ## pseudo-observations
    G <- n

    pcop <- pcopula(cop,g)

    ## compute the test statistic
    stat <- .C(cramer_vonMises_2,
               as.integer(p),
               as.double(u),
               as.integer(n),
               as.double(g),
               as.integer(G),
               as.double(pcop),
               stat = double(1))$stat

    x0 <- u # rcopula(cop,m)

    v <- if (M > 0) rcopula(cop,M) else u

    s0 <- .C(multiplier,
             as.integer(p),
             as.double(x0),
             as.integer(m),
             as.double(g),
             as.integer(G),
             as.double(derCdfWrtParams(cop,g) %*% influCoef(cop,x0,v)),
             as.integer(N),
             s0 = double(N))$s0

    return(list(statistic=stat, pvalue=(sum(s0 >= stat)+0.5)/(N+1),
                parameters=cop@parameters))
}


