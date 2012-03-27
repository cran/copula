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


## CAUTION: sign.fafac is not correct [see below for a corrected version]

## Utility for  polyG()  in ./R/aux-acopula.R :
sign.fafac <- function(alpha, d) { ## also  on (d ; k:= 1:d ; k1:= 0:(d-1))
    stopifnot(length(alpha)==1, 0 < alpha, alpha < 1,
              length(d) == 1, d == round(d), d >= 1)
    k <- 1:d
    k1 <- k - 1L
    s <- unlist(lapply(alpha*k, function(z) prod(z-k1)))
    ss <- sign(s)
    sn <- (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
    ##    --------------------------------------
    stopifnot(ss[s != 0] == sn[s != 0])
    ss
}

sign.fafac(0.8, 4)
sign.fafac(0.8, 5)
sign.fafac(0.8, 6)
sign.fafac(0.8, 16)
sign.fafac(0.8, 17)
## determine signs of the falling factorials
## s <- unlist(lapply(alpha*k, function(z) prod(z-(0:(d-1)))))
## signs <- (-1)^(d-k) * sign(s) ## see  ../misc/sign-polyG.R
## signs  <- (-1)^(d-k) * (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
## signs  =  (-1)^k              * (2*(floor(alpha*k) %% 2) - 1)


## CAUTION: sign.fafac is not correct [see below for a corrected version]
## Ok: It is "proven"
##
##    sign(s) == (-1)^d * (2*(floor(alpha*k) %% 2) - 1)
##
## at whenever sign(s) != 0  {and for the s == 0  case, sign does not matter}

## test 1
alpha <- 0.9
d <- 2
j <- 1:d
signFF(alpha, 1:d, )
sign(choose(alpha*j, d)*(-1)^(d-j))
(-1)^d * (2*(floor(alpha*j) %% 2) - 1) # => formula above is not correct!

## test 2
alpha <- 0.9
d <- 3
j <- 1:d
all(signFF(alpha, 1:d, d) == sign(choose(alpha*j, d)*(-1)^(d-j)))
