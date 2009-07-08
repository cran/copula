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

###################################################################
###################################################################

## EV test based on Cn

evTest <- function(x, r = 6, N = 1000, m = 30, offset=0.5)
{
  nr <- length(r)
  
  ## make pseudo-observations
  p <- ncol(x)
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  if (m > 0)
    {
      y <- seq(1/m, 1 - 1/m, len = m)
      v <- vector("list", p)
      for (i in 1:p)
        v[[i]] <- y
      g <- as.matrix(expand.grid(v))
      m <- nrow(g)
    }
  else 
    {
      g <- u ## pseudo-observations
      m <- n
    }

  ## compute the test statistic
  s <- .C("evtest_stat",
          as.double(u),
          as.integer(n),
          as.integer(p),
          as.double(g),
          as.integer(m),
          as.double(r),
          as.integer(nr),
          as.double(offset),
          stat = double(nr),
          PACKAGE="copula")$stat
  
  s0 <- .C("evtest",
           as.double(u),
           as.integer(n),
           as.integer(p),
           as.double(g),
           as.integer(m),
           as.integer(N),
           as.double(r),
           as.integer(nr),
           as.double(offset),
           s0 = double(N * nr),
           PACKAGE="copula")$s0

  ## combined p-values?
  s0 <- matrix(s0, ncol = nr, byrow = TRUE)
  pval <- apply(s0 >= matrix(s, nrow=N, ncol=nr, byrow=TRUE),
                2, function(x) (sum(x) + 0.5) / (N + 1) )
  comb.s <- sum(s) / nr
  comb.pval <- ( sum( apply(s0, 1, sum) / nr >= comb.s) + 0.5 ) / (N + 1) 

  return(list(statistic=c(s, comb.s),
              pvalue=c(pval, comb.pval), s0=s0))
}

###################################################################
###################################################################
## ev test based on A

evTestA <- function(x, N = 1000, estimator = "CFG", m=30, offset=0.5)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  if (m > 0)
    {
      xis <- yis <- seq(1/m, 1 - 1/m, len = m)
      g <- as.matrix(expand.grid(xis, yis))
      m <- m^2
    }
  else
    {
      g <- u
      m <- n
    }

  ## compute the test statistic
  s <- .C("evtestA_stat",
          as.double(u[,1]),
          as.double(u[,2]),
          as.integer(n),
          as.double(g[,1]),
          as.double(g[,2]),
          as.integer(m),
          as.double(offset),
          as.integer(estimator == "CFG"),
          stat = double(1),
          PACKAGE="copula")$stat
  
  s0 <- .C("evtestA",
           as.double(u[,1]),
           as.double(u[,2]),
           as.integer(n),
           as.double(g[,1]),
           as.double(g[,2]),
           as.integer(m),
           as.double(offset),
           as.integer(estimator == "CFG"),
           as.integer(N),
           s0 = double(N),
           PACKAGE="copula")$s0

  return(list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1),s0=s0))
  
}

