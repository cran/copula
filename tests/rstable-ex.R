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


require(copula)

##  See  limit  alpha --> 1  for fixed gamma :
set.seed(101)
X <- rstable1(1e4, alpha=.9999, beta=1, gamma= .25, delta=1)
summary(X)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##   1592    1592    1593    1594    1593    2092
stopifnot(all(X > 1500))
X <- rstable1(1e4, alpha=.999999, beta=1, gamma= .25, delta=1)
stopifnot(150000 < X, X < 170000)

r1R <- copula:::rstable1R
r1C <- copula:::rstable1C
set.seed(123)

## speed is *very* similar (!):
system.time(Z  <- r1R(10000, .2, beta=1, gamma= 45))
system.time(Z. <- r1C(10000, .2, beta=1, gamma= 45))

ks.test(Z, Z.) # p-value ~ 0.50  they "are the same"

if(!dev.interactive(orNone=TRUE)) pdf("rstable-ex.pdf")
qqplot(Z, Z., log = "xy") # looks nice
acf(log(Z))  # "nice"
acf(log(Z.)) # ditto

cat('Time elapsed: ', proc.time(),'\n') # for ``statistical reasons''
