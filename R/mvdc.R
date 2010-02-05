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


mvdc <- function(copula, margins, paramMargins, marginsIdentical = FALSE) {
  if (marginsIdentical) {
    if (length(margins) == 1) margins <- rep(margins, copula@dimension)
    if (length(paramMargins) == 1) paramMargins <- rep(paramMargins, copula@dimension)
  }
  new("mvdc", copula = copula, margins = margins, paramMargins = paramMargins,
       marginsIdentical = marginsIdentical)
}


## Functions asCall and P0 were kindly supplied by 
## Martin Maechler <maechler@stat.math.ethz.ch>,
## motivated by an application of nor1mix and copula
## from Lei Liu <liulei@virginia.edu>.
## They fixes the function getExpr in the old
## version, which assumed that the parameters to
## [rdpq]<distrib> were vectors.

asCall <- function(fun, param)
{
    cc <-
	if (length(param) == 0)
	    quote(FUN(x))
	else if(is.list(param)) {
	    as.call(c(quote(FUN), c(quote(x), as.expression(param))))
	} else { ## assume that [dpq]<distrib>(x, param) will work
	    as.call(c(quote(FUN), c(quote(x), substitute(param))))
	}
    cc[[1]] <- as.name(fun)
    cc
}

P0 <- function(...) paste(..., sep="")

dmvdc <- function(mvdc, x) {
  dim <- mvdc@copula@dimension
  densmarg <- 1
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    cdf.expr <- asCall(P0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    pdf.expr <- asCall(P0("d", mvdc@margins[i]), mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
    densmarg <- densmarg * eval(pdf.expr, list(x = x[,i]))
  }
  dcopula(mvdc@copula, u) * densmarg
}


pmvdc <- function(mvdc, x) {
  dim <- mvdc@copula@dimension
  if (is.vector(x)) x <- matrix(x, nrow = 1)
  u <- x
  for (i in 1:dim) {
    cdf.expr <- asCall(P0("p", mvdc@margins[i]), mvdc@paramMargins[[i]])
    u[,i] <- eval(cdf.expr, list(x = x[,i]))
  }
  pcopula(mvdc@copula, u)
}


rmvdc <- function(mvdc, n) {
  dim <- mvdc@copula@dimension
  u <- rcopula(mvdc@copula, n)
  x <- u
  for (i in 1:dim) {
    qdf.expr <- asCall(P0("q", mvdc@margins[i]), mvdc@paramMargins[[i]])
    x[,i] <- eval(qdf.expr, list(x = u[,i]))
  }
  x
}
