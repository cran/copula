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


source(system.file("Rsource", "AC-Liouville.R", package="copula"))

require(copula)

### Examples ###################################################################

set.seed(1)

### Archimedean-Simplex copulas ################################################

n <- 1000
theta <- 0.59
d <- 3

U <- rACsimplex(n, d=d, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

pairs(U, gap=0)

hist(U[,1])
hist(U[,2])
hist(U[,3])


### Liouville copulas ##########################################################

## see Figure 3 in McNeil, Neslehova  (2010)

n <- 2000
theta <- 0.6
alpha <- c(1, 5, 20)

U <- rLiouville(n, alpha=alpha, theta=theta, Rdist="Gamma")
cor(U, method="kendall")

pairs(U, gap=0, cex=0.5)

hist(U[,1])
hist(U[,2])
hist(U[,3])


### Archimedean-Liouville copulas ##############################################

## see Figure 4 in McNeil, Neslehova  (2010)

n <- 1000
theta <- 0.59
alpha <- c(1, 3, 4)

U <- rACLiouville(n, alpha=alpha, theta=theta, family="Clayton")
cor(U, method="kendall")

pairs(U, gap=0)

hist(U[,1])
hist(U[,2])
hist(U[,3])

