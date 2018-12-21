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
set.seed(47)

## bivariate comparisons
d <- 2
u <- pobs(matrix(runif(d * m), m, d))

## (Warnings suppressed now via default may.warn=FALSE)
cDer <- rbind(
    clayton = comparederiv(claytonCopula (iTau(claytonCopula(), tau)), u),
    gumbel  = comparederiv(gumbelCopula  (iTau(gumbelCopula(),  tau)), u),
    frank   = comparederiv(frankCopula   (iTau(frankCopula(),   tau)), u),
    plackett= comparederiv(plackettCopula(iTau(plackettCopula(),tau)), u),
    normal  = comparederiv(normalCopula  (iTau(normalCopula(),  tau)), u),
    tC.fixed= comparederiv(tCopula       (iTau(tCopula(), tau), df.fixed = TRUE), u))
cDer
stopifnot(cDer[,"dCdu"      ] <= 0.004, # max: normal   = 0.002166
          cDer[,"dCdtheta"  ] <= 11e-14,# max: tC.fixed = 5.537e-14
          cDer[,"dlogcdu"   ] <= 15e-8, # max: normal   = 7.51e-8
          cDer[,"dlogcdtheta"]<= 6e-9)  # max: normal   = 2.92e-9
showProc.time()


if (doExtras)
{
    ## d-dimensional
    d <- 4 ; set.seed(44)
    u <- pobs(matrix(runif(d * m), m, d))

    nC4 <- normalCopula(rep(iTau(normalCopula(), tau), d * (d-1)/2), dim=d, dispstr = "un")
    tC4 <- tCopula     (rep(iTau(tCopula(),      tau), d * (d-1)/2), dim=d, dispstr = "un",
                        df.fixed = TRUE)
    cD <- rbind(comparederiv(nC4, u),
                comparederiv(tC4, u))
    print(cD, digits = 5)
    stopifnot(apply(cD, 2, max) < c(0.42, 0.18, 2.1e-07, 1.6e-08))
    showProc.time()
}

