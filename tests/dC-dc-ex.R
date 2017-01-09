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
source(system.file("Rsource", "utils.R",     package="copula", mustWork=TRUE))
##-> assertError(), assert.EQ(), ... showProc.time()  +  comparederiv()
showProc.time()

(doExtras <- copula:::doExtras())


m <- 10 # number of random points
tau <- 0.5

## bivariate comparisons
d <- 2
u <- pobs(matrix(runif(d * m), m, d))

## each with  4  warnings  "... numerical differentiation used":
comparederiv(claytonCopula (iTau(claytonCopula(), tau)), u)
comparederiv(gumbelCopula  (iTau(gumbelCopula(),  tau)), u)
comparederiv(frankCopula   (iTau(frankCopula(),   tau)), u)
comparederiv(plackettCopula(iTau(plackettCopula(),tau)), u)
comparederiv(normalCopula  (iTau(normalCopula(),  tau)), u)
comparederiv(tCopula(iTau(tCopula(), tau), df.fixed = TRUE), u)

showProc.time()

if (doExtras)
{
    ## d-dimensional
    d <- 4
    u <- pobs(matrix(runif(d * m), m, d))

    nC4 <- normalCopula(rep(iTau(normalCopula(), tau), d * (d-1)/2), dim=d, dispstr = "un")
    comparederiv(nC4, u)

    tC4 <- tCopula(rep(iTau(tCopula(), tau), d * (d-1)/2), dim=d, dispstr = "un", df.fixed = TRUE)
    comparederiv(tC4, u)
    showProc.time()
}

