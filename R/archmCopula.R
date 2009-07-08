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

archmCopula <- function(family, param, dim = 2, ...) {
  familiesImplemented <- c("clayton", "frank", "gumbel", "amh")
  fam <- pmatch(family, familiesImplemented, -1)
  if (fam == -1)
    stop(paste("Valid family names are", familiesImplemented))
  copula <- switch(fam,
                   claytonCopula(param, dim = dim),
                   frankCopula(param, dim = dim),
                   gumbelCopula(param, dim = dim),
                   amhCopula(param, dim = 2)
                   )
  copula
}



kendallsTauArchmCopula <- function(copula) {
  integrand <- function(x) genFun(copula, x) / genFunDer1(copula, x)
  1 + 4 * integrate(integrand, 0, 1)$value
}

setMethod("kendallsTau", signature("archmCopula"), kendallsTauArchmCopula)
