#################################################################################
##
##   R package Copula by Jun Yan and Ivan Kojadinovic Copyright (C) 2009
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

#################################################################################
## Goodness-of-fit test for extreme value copulas
#################################################################################

## copula is a copula of the desired family whose parameters, if necessary,
## will be used as starting values in fitCopula

gofEVCopula <- function(copula, x, N = 1000, method = "mpl",
                        estimator = "CFG", m = 1000, print.every = 100,
                        optim.method = "Nelder-Mead")
  {
    n <- nrow(x)
    
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    fcop <- fitCopula(copula, u, method, estimate.variance=FALSE,
                      optim.method=optim.method)@copula

    ## where to compute Afun
    g <- seq(0,1-1/m,by=1/m)

    ## compute the test statistic
    s <- .C("cramer_vonMises_Afun",
            as.integer(n),
            as.integer(m),
            as.double(-log(u[,1])),
            as.double(-log(u[,2])),
            as.double(Afun(fcop,g)),
            stat = double(2),
            as.integer(estimator == "CFG"),
            PACKAGE="copula")$stat
    
    ## simulation of the null distribution
    s0 <- matrix(NA, N, 2)
    if (print.every > 0)
      cat(paste("Progress will be displayed every", print.every, "iterations.\n"))
    for (i in 1:N)
      {
        if (print.every > 0 & i %% print.every == 0)
          cat(paste("Iteration",i,"\n"))
        u0 <- apply(rcopula(fcop,n),2,rank)/(n+1)
        
        ## fit the copula
        fcop0 <-  fitCopula(copula, u0, method, estimate.variance=FALSE,
                            optim.method=optim.method)@copula
       
        s0[i,] <-  .C("cramer_vonMises_Afun",
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(Afun(fcop0,g)),
                     stat = double(2),
                     as.integer(estimator == "CFG"),
                     PACKAGE="copula")$stat
      }

    ## corrected version only
    gof <- list(statistic=s[1],
                pvalue=(sum(s0[,1] >= s[1])+0.5)/(N+1),
                parameters=fcop@parameters)

    class(gof) <- "gofCopula"
    gof
  }

#######################################################3
## version for simulations
## was named gofEVCopula before
## not exported

gofAfun <- function(copula, x, N = 1000, method = "mpl", # estimator = "CFG",
                        m = 1000, print.every = 100, optim.method = "Nelder-Mead")
  {
    n <- nrow(x)
    p <- ncol(x)
    
    ## make pseudo-observations
    u <- apply(x,2,rank)/(n+1)

    ## fit the copula
    fcop <- fitCopula(copula, u, method, estimate.variance=FALSE,
                      optim.method=optim.method)@copula

    ## statistic based on Cn
    sCn <- .C("cramer_vonMises",
              as.integer(n),
              as.integer(p),
              as.double(u),
              as.double(pcopula(fcop,u)),
              stat = double(1),
              PACKAGE="copula")$stat
    
    ## where to compute Afun
    g <- seq(0,1-1/m,by=1/m)

    ## compute the CFG test statistic
    sCFG <- .C("cramer_vonMises_Afun",
               as.integer(n),
               as.integer(m),
               as.double(-log(u[,1])),
               as.double(-log(u[,2])),
               as.double(Afun(fcop,g)),
               stat = double(2),
               as.integer(1), # estimator == "CFG"),
               PACKAGE="copula")$stat
    ## compute the Pickands test statistic
    sPck <- .C("cramer_vonMises_Afun",
               as.integer(n),
               as.integer(m),
               as.double(-log(u[,1])),
               as.double(-log(u[,2])),
               as.double(Afun(fcop,g)),
               stat = double(2),
               as.integer(0), # estimator == "CFG"),
               PACKAGE="copula")$stat

    s <- c(sCn, sCFG, sPck)
    
    ## simulation of the null distribution
    s0 <- matrix(NA, N, 5)
    if (print.every > 0)
      cat(paste("Progress will be displayed every", print.every, "iterations.\n"))

    ## set starting values for fitCopula
    copula@parameters <- fcop@parameters
    for (i in 1:N)
      {
        if (print.every > 0 & i %% print.every == 0)
          cat(paste("Iteration",i,"\n"))
        u0 <- apply(rcopula(fcop,n),2,rank)/(n+1)
        
        ## fit the copula
        fcop0 <-  fitCopula(copula, u0, method, estimate.variance=FALSE,
                            optim.method=optim.method)@copula

        sCn0 <- .C("cramer_vonMises",
                   as.integer(n),
                   as.integer(p),
                   as.double(u0),
                   as.double(pcopula(fcop0,u0)),
                   stat = double(1),
                   PACKAGE="copula")$stat

        sCFG0 <-  .C("cramer_vonMises_Afun",
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(Afun(fcop0,g)),
                     stat = double(2),
                     as.integer(1), # estimator == "CFG"),
                     PACKAGE="copula")$stat
        sPck0 <-  .C("cramer_vonMises_Afun",
                     as.integer(n),
                     as.integer(m),
                     as.double(-log(u0[,1])),
                     as.double(-log(u0[,2])),
                     as.double(Afun(fcop0,g)),
                     stat = double(2),
                     as.integer(0), # estimator == "CFG"),
                     PACKAGE="copula")$stat
        s0[i,] <- c(sCn0, sCFG0, sPck0)
      }
    
    return(list(statistic=s,
                pvalue=sapply(1:5, function(i) (sum(s0[,i] >= s[i])+0.5)/(N+1)),
                parameters=fcop@parameters))
  }
