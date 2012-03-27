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


### This function is never used.
### It is commented out on 06/09/2009 because dgamma, which is called by
### sgen, has becomes an intrinsic gfortran function, and gives warnings
### while doing R CMD check with gfortran-4.4.
## rstable <- function(n, alpha, beta, scale = 1, location = 0, iparam = 1) {
##   if (alpha > 2) stop ("alpha must be <= 2")
##   if (beta < -1 | beta > 1) stop("beta must be <= 1 and >= -1")
##   if (scale <= 0) stop("scale must be > 0")
##   val <- double(n)
##   n2 <- n * 2
##   uu <- runif(n2)
##   err <- 0
##   foo <- .Fortran("sgen", as.integer(n), val=as.double(val),
##                   as.double(alpha), as.double(beta),
##                   as.double(scale), as.double(location),
##                   as.integer(n2), as.double(uu),
##                   as.integer(iparam), err=as.integer(err))
##   if (foo$err > 0) stop ("Error in the Frotran code.")
##   foo$val
## }

rPosStable <- function(n, alpha) {
  ## reference: Chambers, Mallows, and Stuck 1976, JASA, p.341
  if (alpha >= 1) stop("alpha must be < 1")
  theta <- runif(n, 0, pi)
  w <- rexp(n)
  a <- sin((1 - alpha) *theta) * (sin(alpha * theta))^(alpha / (1 - alpha)) / (sin(theta))^(1/(1 - alpha))
  (a / w)^((1 - alpha)/alpha)
}
