## Copyright (C) 2016 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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

### Finite Mixtures of Copulas:  mixCopula()
### ========================================

library(copula)
isExplicit <- copula:::isExplicit
(doExtras <- copula:::doExtras())

is32 <- .Machine$sizeof.pointer == 4 ## <- should work for Linux/MacOS/Windows
isMac <- Sys.info()[["sysname"]] == "Darwin"
isSun <- Sys.info()[["sysname"]] == "SunOS"

mC <- mixCopula(list(gumbelCopula(2.5, dim=3),
                     claytonCopula(pi, dim=3),
                     tCopula(0.7, dim=3)),
                c(2,2,4)/8)
mC
stopifnot(dim(mC) == 3, inherits(mC, "mixCopula"))

## mix a Gumbel with a rotated Gumbel (with equal weights 1/2):
mGG <- mixCopula(list(gumbelCopula(2), rotCopula(gumbelCopula(1.5))))
stopifnot(dim(mGG) == 2, inherits(mGG, "mixCopula"),
    all.equal(   rho(mGG), 0.57886340158, tol=1e-10) ,
    all.equal(lambda(mGG), c(lower = 0.206299474016,
                             upper = 0.292893218813), tol = 1e-10)
)
mGG # 4 parameters

RNGversion("3.5.0") # for sample() -- "biased" --> warning ---> *remove* in future!
set.seed(17)
uM  <- rCopula( 600, mC)
uGG <- rCopula(1000, mGG)
## Check  dCopula(*, log= TRUE) ## --- this was __wrong__ for many days (not on CRAN)
stopifnot(
    length(dCopula(uM[1,,drop=FALSE], mC, log=TRUE)) == 1,# was wrong
    all.equal(log(dCopula(uM,  mC)),
		  dCopula(uM,  mC,  log=TRUE), tol = 1e-12),
    all.equal(log(dCopula(uGG, mGG)),
		  dCopula(uGG, mGG, log=TRUE), tol = 1e-12)
)
(llGG1 <- loglikCopula(c(2.5, 1., w = c(2,4)/6), u=uGG, copula = mGG))
(llMC <-  loglikCopula(c(2.5, pi, rho.1=0.7, df = 4, w = c(2,2,4)/8), u = uM, copula = mC))
## discrepancy 32 bit <-> 64 bit --- still (after dMixCopula() bug fix):
stopifnot(
    all.equal(llGG1, 177.452426, ## 32 bit (Windows, Linux(FC 24)): 188.0358
              tol = if(isSun || isMac || is32) 0.08 else 7e-7),
    all.equal(llMC,  532.8757887, ## 32 bit: 551.8439
              tol = if(isSun || isMac || is32) 0.05 else 7e-7)
)

## "free" weights --- estimation == FIXME: will change after re-parametrization
optCtrl <- list(maxit = 1000, trace = TRUE)
if (doExtras) { # slowish
    st0 <- system.time(
        f. <- fitCopula(mC, uM, optim.method = "BFGS", optim.control=optCtrl, traceOpt=TRUE))
    ## converges, .. no var-cov
    ## (which would fail with Error in  solve.default(Sigma.n, t(dlogcdtheta(copula, u) - S)) :
    ##            system is computationally singular: reciprocal condition number = 3.88557e-17
    print(st0) #
    print(lf. <- logLik(f.))
    print(summary(f.))
}


if (doExtras) { # slowish
    st1 <- system.time(
        ff <- fitCopula(mC, uM, optim.method = "Nelder-Mead", optim.control=optCtrl))
    ## converges (vcov : Error in solve(..)... (rec.cond.number 4.783e-17)
    print(st1) # 11 sec
    print(lff <- logLik(ff))
    print(summary(ff))
}


if (doExtras) { # slowish
    st2 <- system.time(
        f2 <- fitCopula(mC, uM, optim.method = "L-BFGS-B",    optim.control=optCtrl))
    ## converges (vcov: Error in solve(..)... (rec.cond.number 3.691e-17)
    print(st2) # 28 sec
    print(lf2 <- logLik(f2))
    print(summary(f2))
}

## partially fixed

(tX4 <- tCopula( 0.2, df = 5, df.fixed=TRUE))
(tn3 <- tCopula(-0.5, df = 3))
getTheta(tX4, attr = TRUE) # freeOnly = TRUE is default
## --> *not* showing df=5 as it is fixed
(m3 <- mixCopula(list(normalCopula(0.4), tX4, tn3), w = (1:3)/6))
                                        # -> shows 'm2.df := 5' as fixed !
th. <- getTheta(m3, attr = TRUE)# ditto
(th <- getTheta(m3, named= TRUE))
trueAt <- function(i, n) { r <- logical(n); r[i] <- TRUE; r }
stopifnot(
    identical(th, structure(th., param.lowbnd=NULL, param.upbnd=NULL)),
    identical(names(th), c("m1.rho.1", "m2.rho.1", "m3.rho.1", "m3.df",
                           paste0("w", 1:3)))
    ,
    identical(isFree(m3), !trueAt(3, n=8))# free every but at [3]
)

fixedParam(m3) <- trueAt(c(3, 8), n=8)
