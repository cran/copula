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


tC2.F <- tCopula(df=2, df.fixed=TRUE)
(cm <- setTheta(tC2.F, value=-0.5, freeOnly=TRUE)) ## setTheta() had failed
(cp <- setTheta(tC2.F, value= 0.5, freeOnly=TRUE))
(c2 <- setTheta(tC2.F, value=c(0.5,3), freeOnly=FALSE)) ## had failed
stopifnot(all.equal(getTheta(cm), -0.5), all.equal(getTheta(cp), +0.5),
          all.equal(getTheta(cm, freeOnly=FALSE, named=TRUE),
                    c(rho.1 = -0.5, df = 2)),
          all.equal(getTheta(cp, freeOnly=FALSE), c(0.5, 2)),
          all.equal(getTheta(c2, freeOnly=FALSE, named=TRUE),
                    c(rho.1 = 0.5, df = 3)))

(N3 <- normalCopula(c(0.5,0.3,0.2), dim=3, dispstr = "un"))
fixedParam(N3) <- c(TRUE, FALSE, FALSE); N3
(N3.2 <- setTheta(N3,     c( 0.4, 0.2) -> t2)) # partially fixed: works since 2017-02-11
(N3.3 <- setTheta(N3, c(0.6, 0.4, 0.2) -> t3, freeOnly=FALSE))
stopifnot(all.equal(getTheta(N3.2), t2),
	  all.equal(getTheta(N3.3), t3[-1]),
	  all.equal(getTheta(N3.2, freeOnly=FALSE), c(0.5, t2)),
	  all.equal(getTheta(N3.3, freeOnly=FALSE), t3))

(tC3 <- tCopula(c(0.5,0.3,0.2), dim=3, dispstr = "un")) # df = 4 is not fixed
fixedParam(tC3) <- c(TRUE, FALSE, FALSE, TRUE); tC3 #-> df fixed, too
(tC3.2 <- setTheta(tC3,     c( 0.4, 0.1)    -> t2)) # partially fixed: works since 2017-02-11
(tC3.3 <- setTheta(tC3, c(0.6, 0.4, 0.1, 3) -> t3, freeOnly=FALSE))
stopifnot(all.equal(getTheta(tC3.2), t2),
	  all.equal(getTheta(tC3.3), t3[-c(1,4)]),
	  all.equal(getTheta(tC3.2, freeOnly=FALSE), c(0.5, t2, 4)),
	  all.equal(getTheta(tC3.3, freeOnly=FALSE), t3))

fixedParam(tC3) <- c(TRUE, FALSE, FALSE, FALSE); tC3 # df remains free
(tC3u.2 <- setTheta(tC3,     c( 0.4, 0.2, 5) -> t2)) # partially fixed: works since 2017-02-11
(tC3u.3 <- setTheta(tC3, c(0.6, 0.4, 0.2, 3) -> t3, freeOnly=FALSE))
stopifnot(all.equal(getTheta(tC3u.2), t2),
	  all.equal(getTheta(tC3u.3), t3[-1]),
	  all.equal(getTheta(tC3u.2, freeOnly=FALSE), c(0.5, t2)),
	  all.equal(getTheta(tC3u.3, freeOnly=FALSE), t3))



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

    comparederiv(nc2, v) ## from  ../inst/Rsource/utils.R
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

