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


archmCopula <- function(family, param, dim = 2L, ...) {
  familiesImplemented <- c("clayton", "frank", "gumbel", "amh")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  dim <- as.integer(dim)
  switch(fam,
         claytonCopula(param, dim = dim),
         frankCopula  (param, dim = dim),
         gumbelCopula (param, dim = dim),
         amhCopula    (param, dim = 2L)
         )
}


kendallsTauArchmCopula <- function(copula) {
  integrand <- function(x) genFun(copula, x) / genFunDer1(copula, x)
  1 + 4 * integrate(integrand, 0, 1)$value
}

setMethod("kendallsTau", signature("archmCopula"), kendallsTauArchmCopula)
