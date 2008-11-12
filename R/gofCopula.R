#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2008
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

## Goodness-of-fit test for copulas: wrapper function
## Calls either gofMCLT.PL, gofMCLT.KS or gofPB
## See also the file mgofMclt.R

## copula is a copula of the desired family

gofCopula <- function(copula, x, N = 1000, method = "mpl",
                      simulation = "pb", grid = "h0", R = 10,
                      m = nrow(x), G = nrow(x), M = 2500)
  {
    x <- as.matrix(x)
    
    n <- nrow(x)
    p <- ncol(x)
    
    if (p < 2)
      stop("The data should be at least of dimension 2")
    
    if (n < 2)
      stop("There should be at least 2 observations")
    
    if (copula@dimension != p)
      stop("The copula and the data should be of the same dimension")


    if (simulation == "pb") ## parametric bootstrap
      gofPB(copula, x, N, method)
    else if (simulation == "mult") ## multiplier
      {
        if (method == "mpl")
          gofMCLT.PL(copula, x, N, m, grid, G, R, M)
        else if (method %in% c("irho","itau"))
          {
            if (copula@dimension != 2)
              stop("The simulation method 'mult' can be used in combination with the estimation methods 'irho' and 'itau' only in the bivariate case.") 
            gofMCLT.KS(copula, x, N, m, method, grid, G, R, M)
          }
        else stop("Invalid estimation method")
      }
    else stop("Invalid simulation method")
  }

##############################################################################

## Goodness-of-fit test based on the parametric bootstrap
## as proposed by Genest et al. (2008)

## copula is a copula of the desired family whose parameters, if necessary,
## will be used as starting values in fitCopula

gofPB <- function(copula, x, N, method, P=100)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    fcop <- fitCopula(u,copula,method,estimate.variance=FALSE)@copula

    ## compute the test statistic
    s <- .C("cramer_vonMises",
            as.integer(n),
            as.integer(p),
            as.double(u),
            as.double(pcopula(fcop,u)),
            stat = double(1),
            PACKAGE="copula")$stat
    
    ## simulation of the null distribution
    s0 <- numeric(N)
    cat(paste("Progress will be displayed every",P,"iterations.\n"))
    for (i in 1:N)
      {
        if (i %% P == 0)
          cat(paste("Iteration",i,"\n"))
        u0 <- apply(rcopula(fcop,n),2,rank)/(n+1)
        
        ## fit the copula
        fcop0 <-  fitCopula(u0,copula,method,estimate.variance=FALSE)@copula
        
        s0[i] <- .C("cramer_vonMises",
                    as.integer(n),
                    as.integer(p),
                    as.double(u0),
                    as.double(pcopula(fcop0,u0)),
                    stat = double(1),
                    PACKAGE="copula")$stat
      }
    
    return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1),
                parameters=fcop@parameters))
  }

##############################################################################
## additional influence terms

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

##############################################################################

## Goodness-of-fit test based on the multiplier approach
## and rank correlation coefficients

gofMCLT.KS <- function(cop, x, N, m, method, grid, G, R, M)
  {
    n <- nrow(x)
    p <- 2
  
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    cop <- fitCopula(u,cop,method,estimate.variance=FALSE)@copula
    
    ## perform R replications
    stat <- pval <- numeric(R)
    for (i in 1:R)
      { 
        
        ## grid points where to evaluate the process
        if (grid == "h0")  ## from the H0 copula 
          g <- rcopula(cop,G)
        else if (grid == "po") ## pseudo-observations
          {
            g <- u
            G <- n
          }
        else stop("Invalid grid")
        
        pcop <- pcopula(cop,g)
        
        ## compute the test statistic
        stat[i] <- .C("cramer_vonMises_2",
                      as.integer(p),
                      as.double(u),
                      as.integer(n),
                      as.double(g),
                      as.integer(G),
                      as.double(pcop),
                      stat = double(1),
                      PACKAGE="copula")$stat
        
        ## generate realizations under H0
        x0 <- rcopula(cop,m) ## default
        
        ## prepare influence coefficients              
        if (method == "itau") ## kendall's tau
          influ <- 4 * (2 * pcopula(cop,x0) - x0[,1] - x0[,2] + (1 - kendallsTau(cop))/2) / tauDer(cop)
        else if (method == "irho") ## Spearman's rho
          {
            ## integrals computed from M realizations by Monte Carlo
            y0 <- rcopula(cop,M)
            influ <- (12 * (x0[,1] * x0[,2] + influ.add(x0, y0, y0[,2],y0[,1])) -
                      3 - spearmansRho(cop)) / rhoDer(cop)
          }
        
        ## Simulate under H0
        s0 <- .C("multiplier",
                 as.integer(p),
                 as.double(x0),
                 as.integer(m),
                 as.double(g), 
                 as.integer(G), 
                 as.double(pcop), 
                 as.double(derCdfWrtArgs(cop,g)),
                 as.double(derCdfWrtParams(cop,g) %*% influ),
                 as.integer(N),
                 s0 = double(N),
                 PACKAGE="copula")$s0
    
        pval[i] <- (sum(s0 >= stat[i])+0.5)/(N+1)
        
      }

    return(list(statistic=median(stat), pvalue=median(pval),
                sd.pvalues=sd(pval), parameters=cop@parameters))
  }
