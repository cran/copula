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


### Computing probabilities of falling in hyperrectangles

### prob() --- Generic and all Methods here
### ======

##' @title Compute the probability P[l < U <= u]  where U ~ copula x
##' @param x copula object
##' @param l d-vector of lower "integration" limits
##' @param u d-vector of upper "integration" limits
##' @return the probability that a random vector following the given copula
##'         falls in the hypercube with lower and upper corner l and u, respectively.
##' @author Marius Hofert, Martin Maechler
setGeneric("prob", function(x, l, u) standardGeneric("prob"))

setMethod("prob", signature(x ="Copula"),
          function(x, l,u) {
              d <- dim(x)
              stopifnot(is.numeric(l), is.numeric(u),
                        length(u) == d, d == length(l),
                        0 <= l, l <= u, u <= 1)
              if(d > 30)
		  stop("prob() for copula dimensions > 30 are not supported (yet)")
              D <- 2^d
              m <- 0:(D - 1)
              ## digitsBase() from package 'sfsmisc' {slightly simplified} :
              ## Purpose: Use binary representation of 0:N
              ## Author: Martin Maechler, Date:  Wed Dec  4 14:10:27 1991
              II <- matrix(0, nrow = D, ncol = d)
              for (i in d:1L) {
                  II[,i] <- m %% 2L + 1L
                  if (i > 1) m <- m %/% 2L
              }
              ## Sign: the ("u","u",...,"u") case has +1; = c(2,2,...,2)
              Sign <- c(1,-1)[1L + (- rowSums(II)) %% 2]
              U <- array(cbind(l,u)[cbind(c(col(II)), c(II))], dim = dim(II))
              sum(Sign * pCopula(U, x))
          })
