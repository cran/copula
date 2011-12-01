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
# Test of exchangeability based on An
###################################################################

exchEVTest <- function(x, N = 1000,  estimator = "CFG", derivatives = "Cn", m = 100)
{
  ## make pseudo-observations
  n <- nrow(x)
  u <- apply(x,2,rank)/(n+1)

  ## make grid
  g <- seq(1/m, 0.5, len = m)

  ## compute the test statistic
  s <- .C("evsymtest_stat",
          as.double(-log(u[,1])),
          as.double(-log(u[,2])),
          as.integer(n),
          as.double(g),
          as.integer(m),
          as.integer(estimator == "CFG"),
          stat = double(1),
          PACKAGE="copula")$stat
  
  if (derivatives == "Cn")
    s0 <- .C("evsymtest",
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N),
             PACKAGE="copula")$s0
  else
    s0 <- .C("evsymtest_derA",
             as.double(u[,1]),
             as.double(u[,2]),
             as.integer(n),
             as.double(g),
             as.integer(m),
             as.integer(estimator == "CFG"),
             as.integer(N),
             s0 = double(N),
             PACKAGE="copula")$s0
  
  excht <- list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1))
  class(excht) <- "exchTest"
  excht
}

print.exchTest <- function(x, ...)
{
  cat("Statistic:", x$statistic,
      "with p-value", x$pvalue, "\n\n")
}

###################################################################
## Test of exchangeability based on Cn
###################################################################

exchTest <- function(x, N = 1000, m = 0)
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
  s <- .C("exchtestCn_stat",
          as.double(u[,1]),
          as.double(u[,2]),
          as.integer(n),
          as.double(g[,1]),
          as.double(g[,2]),
          as.integer(m),
          stat = double(1),
          PACKAGE="copula")$stat
  
  s0 <- .C("exchtestCn",
           as.double(u[,1]),
           as.double(u[,2]),
           as.integer(n),
           as.double(g[,1]),
           as.double(g[,2]),
           as.integer(m),
           as.integer(N),
           s0 = double(N),
           PACKAGE="copula")$s0

  excht <- list(statistic=s, pvalue=(sum(s0 >= s)+0.5)/(N+1))
  class(excht) <- "exchTest"
  excht
}
