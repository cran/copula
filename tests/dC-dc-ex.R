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

(doExtras <- copula:::doExtras())

if (doExtras)
{
    m <- 10 # number of random points
    tau <- 0.5

    ## Returns max error
    comparederiv <- function(cop, u) {

        c(dCdu = max(abs((copula:::dCdu(cop, u) -
                          copula:::dCduCopulaNum(cop, u)))),
          dCdtheta = max(abs(copula:::dCdtheta(cop, u) -
                             copula:::dCdthetaCopulaNum(cop, u))),
          dlogcdu = max(abs(copula:::dlogcdu(cop, u) -
                            copula:::dlogcduCopulaNum(cop, u))),
          dlogcdtheta = max(abs(copula:::dlogcdtheta(cop, u) -
                                copula:::dlogcdthetaCopulaNum(cop, u))))
    }



    ## bivariate comparisons
    d <- 2
    u <- pobs(matrix(runif(d * m), m, d))

    cop  <- claytonCopula(iTau(claytonCopula(), tau))
    comparederiv(cop, u)

    cop  <- gumbelCopula(iTau(gumbelCopula(), tau))
    comparederiv(cop, u)

    cop  <- frankCopula(iTau(frankCopula(), tau))
    comparederiv(cop, u)

    cop  <- plackettCopula(iTau(plackettCopula(), tau))
    comparederiv(cop, u)

    cop  <- normalCopula(iTau(normalCopula(), tau))
    comparederiv(cop, u)

    cop  <- tCopula(iTau(tCopula(), tau), df.fixed = TRUE)
    comparederiv(cop, u)

    ## d-dimensional
    d <- 4
    u <- pobs(matrix(runif(d * m), m, d))

    cop <- normalCopula(rep(iTau(normalCopula(), tau), d * (d-1)/2), dim=d, dispstr = "un")
    comparederiv(cop, u)

    cop <- tCopula(rep(iTau(tCopula(), tau), d * (d-1)/2), dim=d, dispstr = "un", df.fixed = TRUE)
    comparederiv(cop, u)
}

