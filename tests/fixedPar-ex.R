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

showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})

### TEST FITTING ##########################################################

    n <- 100
    ## with normal copulas
    nc3  <- normalCopula(dim = 3, c(.6,.3,.2), dispstr = "un")
    nc3@parameters
    set.seed(4521)
    x <- rCopula(n, nc3)
    u <- pobs(x)

    fitCopula(nc3, data = u)
    fitCopula(nc3, data = u, estimate.variance = FALSE)
    fitCopula(nc3, data = x, method = "ml")
    fitCopula(nc3, data = x, method = "ml", estimate.variance = FALSE)
    fitCopula(nc3, data = u, method = "itau")
    fitCopula(nc3, data = u, method = "itau", estimate.variance = FALSE)
    fitCopula(nc3, data = u, method = "irho")
    fitCopula(nc3, data = u, method = "irho", estimate.variance = FALSE)

showProc.time()

    nc2  <- normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                         dispstr = "un")
    nc2@parameters

    fitCopula(nc2, data = u)
    fitCopula(nc2, data = u, estimate.variance = FALSE)
    fitCopula(nc2, data = x, method = "ml")
    fitCopula(nc2, data = x, method = "ml", estimate.variance = FALSE)
    fitCopula(nc2, data = u, method = "itau")
    fitCopula(nc2, data = u, method = "itau", estimate.variance = FALSE)
    fitCopula(nc2, data = u, method = "irho")
    fitCopula(nc2, data = u, method = "irho", estimate.variance = FALSE)

showProc.time()


    nc1  <- normalCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                         dispstr = "un")
    nc1@parameters

    fitCopula(nc1, data = u)
    fitCopula(nc1, data = x, method = "ml")
    fitCopula(nc1, data = u, method = "itau")
    fitCopula(nc1, data = u, method = "irho")

showProc.time()


    ## with t copulas (df.fixed = FALSE)
    tc3df  <- tCopula(dim = 3, c(.6,.3,.2), dispstr = "un")
    tc3df@parameters
    set.seed(4521)
    x <- rCopula(n, tc3df)
    u <- pobs(x)
    fitCopula(tc3df, data = u)
    fitCopula(tc3df, data = x, method = "ml")
    fitCopula(tc3df, data = u, method = "itau")
    fitCopula(tc3df, data = u, estimate.variance = FALSE)
    fitCopula(tc3df, data = x, method = "ml", estimate.variance = FALSE)
    fitCopula(tc3df, data = u, method = "itau", estimate.variance = FALSE)
    fitCopula(tc3df, data = u, method = "itau.mpl")

showProc.time()


    tc2df  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                    dispstr = "un")
    tc2df@parameters

    fitCopula(tc2df, data = u)
    fitCopula(tc2df, data = x, method = "ml")
    fitCopula(tc2df, data = u, method = "itau")
    fitCopula(tc2df, data = u, estimate.variance = FALSE)
    fitCopula(tc2df, data = x, method = "ml", estimate.variance = FALSE)
    fitCopula(tc2df, data = u, method = "itau", estimate.variance = FALSE)
    fitCopula(tc2df, data = u, method = "itau.mpl")
    ## fitCopula(tc2df, data = u, method = "irho")
showProc.time()


    tc1df  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                    dispstr = "un")
    tc1df@parameters

    fitCopula(tc1df, data = u)
    fitCopula(tc1df, data = x, method = "ml")
    fitCopula(tc1df, data = u, method = "itau")
    fitCopula(tc1df, data = u, method = "itau.mpl")
    ## fitCopula(tc1df, data = u, method = "irho")
showProc.time()


    ## with t copulas (df.fixed = TRUE)
    tc2  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, FALSE, FALSE)),
                    dispstr = "un", df.fixed = TRUE)
    tc2@parameters

    fitCopula(tc2, data = u)
    fitCopula(tc2, data = u, method = "itau")
    ## fitCopula(tc2, data = u, method = "irho")
showProc.time()


    tc1  <- tCopula(dim = 3, fixParam(c(.6,.3,.2), c(TRUE, TRUE, FALSE)),
                    dispstr = "un", df.fixed = TRUE)
    tc1@parameters

    fitCopula(tc1, data = u)
    fitCopula(tc1, data = u, method = "itau")
    ##fitCopula(tc1, data = u, method = "irho")
showProc.time()

### TEST dC-dc functions #####################################################

    ## d*du functions should return the same result as when unfixed
    ## d*dtheta functions should return "columns" corresponding to free params
    testdCdc <- function(cop, v, cop.unfixed) {
        fixed <- attr(cop@parameters, "fixed")
        if (.hasSlot(cop, "df.fixed")) fixed <- fixed[-length(fixed)]
        stopifnot(all.equal(copula:::dCdu(cop, v), copula:::dCdu(cop.unfixed, v)),
                  all.equal(copula:::dCdtheta(cop, v),
                            copula:::dCdtheta(cop.unfixed, v)[, !fixed, drop = FALSE]),
                  all.equal(copula:::dlogcdu(cop, v), copula:::dlogcdu(cop.unfixed, v)),
                  all.equal(copula:::dlogcdtheta(cop, v),
                            copula:::dlogcdtheta(cop.unfixed, v)[, !fixed, drop = FALSE]))
    }

    ## random points in unit cube
    set.seed(7615)
    v <- matrix(runif(15), 5, 3)

    ## normal
    testdCdc(nc2, v, nc3)
    testdCdc(nc1, v, nc3)

    ## t with df.fixed = TRUE
    tc3  <- tCopula(dim = 3, c(.6,.3,.2), dispstr = "un", df.fixed = TRUE)
    testdCdc(tc2, v, tc3)
    testdCdc(tc1, v, tc3)

showProc.time()

    ## Compare true and numerical derivatives
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
    comparederiv(nc2, v)
    comparederiv(nc1, v)
    comparederiv(tc2, v)
    comparederiv(tc1, v)

showProc.time()

if (doExtras) {## ~ 7 secs

### Multiplier GOF #####################################################

    ## check size of mult GOF test briefly
    do1 <- function(n, cop) {
        u <- pobs(rCopula(n, cop))
        gofCopula(cop, pobs(u), sim = "mult")$p.value
    }
    n <- 100
    M <- 10 #1000
    mM <- sapply(list(nc2=nc2, nc1=nc1, tc2=tc2, tc1=tc1
                    ## , tc2df, tc1df
                      ), function(COP) mean(replicate(M, do1(n, COP))))
    print(mM)
    print(mM < 0.05)

showProc.time()
}

