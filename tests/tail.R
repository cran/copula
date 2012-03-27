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
## require(fCopulae)

numTailIndexLower <- function(copula, u) {
  ## u is a vector approaching 0
  pcopula(copula, cbind(u, u)) / u
}

numTailIndexUpper <- function(copula, u) {
  #3 u is a vector approaching 1
  (1 - 2 * u + pcopula(copula, cbind(u, u))) / (1 - u)
}

# Upper Tail Dependence
s = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
S = cbind(s, s)

# R/Copula:
C = pcopula(copula = gumbelCopula(param=40, dim = 2), S)
C1 = (1-2*s+C)/(1-s)

# SPlus/Finmetrics
# C = archm.copula(family="gumbel", param=20)
# C = pcopula(C, s, s)
C = c(0.8615672, 0.9300288, 0.9858872, 0.9929363, 0.9985861,
    0.9992930, 0.9998586)
C2 = (1-2*s+C)/(1-s)

## # R/Rmetrics:
## C = parchmCopula(S, alpha=40, type = 4, alternative = TRUE)
## C3 = (1-2*s+C)/(1-s)
## cbind(s, C1, C2, C3)


# Lower Tail Dependence
s = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
s = 1 - s
S = cbind(s, s)

# R/Copula:
C = pcopula(copula = gumbelCopula(param=20, dim = 2), S)
C1 = C/s

# SPlus/Finmetrics:
# C = archm.copula(family="gumbel", param=2)
# pcopula(C, s, s)
C = c(0.03852888470, 0.01445658570, 0.00148447496, 0.00055699612,
    0.00005719516, 2.146044e-005, 2.203666e-006)
C2 = C/s

## # R/Rmetrics:
## C = parchmCopula(S, alpha=20, type = 4, alternative = FALSE)
## C3 = C/s
## cbind(s, C1, C2, C3)


# Upper Tail Dependence
s = c(0.90, 0.95, 0.99, 0.995, 0.999, 0.9995, 0.9999)
S = cbind(s, s)

# R/Copula:
C = pcopula(copula = frankCopula(param=2, dim = 2), S)
C1 = (1-2*s+C)/(1-s)

# SPlus/Finmetrics
# C = archm.copula(family="gumbel", param=2)
# pcopula(C, s, s)
C = c(0.8615672, 0.9300288, 0.9858872, 0.9929363, 0.9985861,
    0.9992930, 0.9998586)
C2 = (1-2*s+C)/(1-s)

## # R/Rmetrics:
## C = parchmCopula(S, alpha=2, type = 5, alternative = TRUE)
## C3 = (1-2*s+C)/(1-s)
## cbind(s, C1, C2, C3)


