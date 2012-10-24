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

##' Generate  stable(alpha, beta=1, pm=1)  random numbers
##'
##' only called for rCopula(<Gumbel>, .)
##' @title Generate  stable(alpha, beta=1, pm=1)  random numbers
##' @param n integer
##' @param alpha number in (0, 1)l
##' @references: Chambers, Mallows, and Stuck 1976, JASA, p.341, formula (2.2)
##' @note: The formula below is (somewhat, for small alpha) *WRONG*:
##'   it lacks a multiplicative factor   cos(pi*alpha/2)^-(1/alpha)
##'   as can be proven by direct comparison with  rstable1() in ./rstable1.R
##' @return numeric vector of length n
rPosStable <- function(n, alpha) {
  if (alpha >= 1) stop("alpha must be < 1")
  theta <- runif(n, 0, pi)
  W <- rexp(n)
  ## a <- sin((1 - alpha) *theta) * sin(alpha * theta)^(alpha / (1 - alpha)) /
  ##     sin(theta)^(1/(1 - alpha))
  a <- sin((1 - alpha) *theta) *
      (sin(alpha * theta)^alpha / sin(theta)) ^ (1/(1 - alpha))
  (a / W)^((1 - alpha)/alpha)
}

