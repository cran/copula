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

## copula is a copula of the desired family whose parameters, if necessary,
## will be used as starting values in fitCopula

gofCopula <- function(copula, x, N=10000, estimation="likelihood", method="parametric.bootstrap", M=2500)
  {
    x <- as.matrix(x)
    
    n <- nrow(x)
    p <- ncol(x)
    
    if (p < 2) stop("The data should be at least of dimension 2")
    
    if (n < 2) stop("There should be at least 2 observations")
    
    if (copula@dimension != p) stop("The copula and the data should be of the same dimension")

    ## how should the parameters be estimated?
    ESTIMATION <- c("likelihood", "kendall", "spearman")
    estimation <-  pmatch(estimation, ESTIMATION)
    if(is.na(estimation))
      stop("Invalid estimation method")
    if(estimation == -1)
      stop("Ambiguous estimation method")

    if ((estimation == 2 || estimation == 3) && p != 2)
      stop("Estimation based on Kendall's tau or Spearman's rho requires p = 2")
    
    if ((estimation == 2 || estimation == 3) && length(copula@parameters) != 1)
      stop("Estimation based on Kendall's tau or Spearman's rho will work only for one-parameter copulas. If you are using the t-copula, use 'df.fixed=TRUE'")

    ## how should the test statistic be simulated under the null?
    METHODS <- c("multiplier", "parametric.bootstrap")
    method <-  pmatch(method, METHODS)
    if(is.na(method))
      stop("Invalid simulation method")
    if(method == -1)
      stop("Ambiguous simulation method")

    ## call the appropriate function
    if (method == 1 && estimation == 1)
      gofMCLT.PL(copula, x, N, M)
    else if (method == 1 && estimation > 1)
      gofMCLT.KS(copula, x, N, estimation, M)
    else if (method == 2)
      gofPB(copula, x, N, estimation)
    else stop("Invalid simulation or estimation method")
    
  }

##############################################################################

## Goodness-of-fit test based on the parametric bootstrap
## as proposed by Genest et al. (2008)

## copula is a copula of the desired family whose parameters, if necessary,
## will be used as starting values in fitCopula

gofPB <- function(copula, x, N, estimation)
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    if (estimation==1)
      fcop <- fitCopula(u,copula,copula@parameters)@copula
    else if (estimation==2)
      fcop <- fitCopulaKendallsTau(copula,u)
    else if (estimation==3)
      fcop <- fitCopulaSpearmansRho(copula,u)

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
    for (i in 1:N)
      {
        cat(paste("Iteration",i,"\n"))
        x0 <- rcopula(fcop,n)
        u0 <- apply(x0,2,rank)/(n+1)
        
         ## fit the copula
        if (estimation==1)
          fcop0 <- fitCopula(u0,copula,fcop@parameters)@copula
        else if (estimation==2)
          fcop0 <- fitCopulaKendallsTau(copula,u0)
        else if (estimation==3)
          fcop0 <- fitCopulaSpearmansRho(copula,u0)
        
        s0[i] <- .C("cramer_vonMises",
                    as.integer(n),
                    as.integer(p),
                    as.double(u0),
                    as.double(pcopula(fcop0,u0)),
                    stat = double(1),
                    PACKAGE="copula")$stat
      }
    
    return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1)))
  }

##############################################################################

fitCopulaKendallsTau <- function(cop,u)
{
    cop@parameters <- calibKendallsTau(cop,cor(u[,1],u[,2],method="kendall"))
    return(cop)
}

##############################################################################

fitCopulaSpearmansRho <- function(cop,u)
{
    cop@parameters <- calibSpearmansRho(cop,cor(u[,1],u[,2],method="spearman"))
    return(cop)
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

gofMCLT.KS <- function(cop, x, N, estimation, M)
  {
    n <- nrow(x)
    p <- 2
  
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    if (estimation==2)
      cop <- fitCopulaKendallsTau(cop,u)
    else if (estimation==3)
      cop <- fitCopulaSpearmansRho(cop,u)
    else stop("Invalid estimation method")
    
    ## estimate of theta
    alpha <-  cop@parameters
    
    ## grid points where to evaluate the process
    g <- rcopula(cop,n) ## default 
    pcop <- pcopula(cop,g)
    
    ## compute the test statistic
    s <- .C("cramer_vonMises_2",
            as.integer(p),
            as.double(u),
            as.integer(n),
            as.double(g),
            as.integer(n),
            as.double(pcop),
            stat = double(1),
            PACKAGE="copula")$stat
   
    ## generate realizations under H0
    x0 <- rcopula(cop,n) ## default
      
    ## prepare influence coefficients              
    if (estimation==2) ## kendall's tau
      influ <- 4 * (2 * pcopula(cop,x0) - x0[,1] - x0[,2] + (1 - kendallsTau(cop))/2) / tauDer(cop)
    else if (estimation==3) ## Spearman's rho
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
             as.integer(n),
             as.double(g), 
             as.integer(n), 
             as.double(pcop), 
             as.double(derCdfWrtArgs(cop,g)),
             as.double(derCdfWrtParams(cop,g) %*% influ),
             as.integer(N),
             s0 = double(N),
             PACKAGE="copula")$s0


    return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1), parameters=alpha))
  }
