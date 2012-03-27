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


x <- rcopula(claytonCopula(3), 200)

## Does the Gumbel family seem to be a good choice?
gofCopula(gumbelCopula(1), x)
## What about the Clayton family?
gofCopula(claytonCopula(1), x)

## The same with a different estimation method
gofCopula(gumbelCopula(1), x, method="itau")
gofCopula(claytonCopula(1), x, method="itau")


## A three-dimensional example
x <- rcopula(tCopula(c(0.5, 0.6, 0.7), dim = 3, dispstr = "un"),200)

## Does the Clayton family seem to be a good choice?
gofCopula(gumbelCopula(1, dim = 3), x)
## What about the t copula?
t.copula <- tCopula(rep(0, 3), dim = 3, dispstr = "un", df.fixed=TRUE)
gofCopula(t.copula, x)

## The same with a different estimation method
gofCopula(gumbelCopula(1, dim = 3), x, method="itau")
gofCopula(t.copula, x, method="itau")

## The same using the multiplier approach
gofCopula(gumbelCopula(1, dim = 3), x, simulation="mult")
gofCopula(t.copula, x, simulation="mult")
