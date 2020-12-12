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
(doExtras <- copula:::doExtras() && getRversion() >= "3.4") # so have withAutoprint(.)

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
mGG # 4 parameters
stopifnot(exprs = {
    dim(mGG) == 2
    inherits(mGG, "mixCopula")
    is(mGG,"mixExplicitCopula") # ==> Jscore() works  ==> var.mpl()
    all.equal(   rho(mGG), 0.57886340158, tol=1e-10)
    all.equal(lambda(mGG), c(lower = 0.206299474016,
                             upper = 0.292893218813), tol = 1e-10)
})


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

## "free" weights -------------

optCtrl <- list(maxit = 1000, trace = TRUE)
## The real proof: using "arbitrary"  initial parameters (not very close to truth):
mC. <- setTheta(mC, c(1, 0.5, rho=0, df=7,   w = c(1,1,1)/3))

if(doExtras) withAutoprint({
    st0 <- system.time(
        f. <- fitCopula(mC, uM, optim.method = "BFGS", optim.control=list(maxit=999), traceOpt=TRUE))
    ## converges, .. no var-cov
    ## (which would fail with Error in  solve.default(Sigma.n, t(dlogcdtheta(copula, u) - S)) :
    ##            system is computationally singular: reciprocal condition number = 3.88557e-17
    st0 # was 4.4, now 9.3 sec (nb-mm4)
    coef(f.) ; logLik(f.)
    coef(f., orig=FALSE) # -> l-scale
    summary(f.)
    stopifnot(exprs = {
        all.equal(logLik(f.), structure(536.93419, nobs = 600L, df = 7L, class = "logLik"),
                  tol = 8e-10) # 7.206e-11)
        })
    ## using 'mC.' which "does not know" the true parameters (for 'start');
    ## now works thanks to new getInitParam() "mixCopula" method :
    system.time(
        f.w <- fitCopula(mC., uM, optim.method = "BFGS", optim.control=list(maxit=999), traceOpt=TRUE))
     ## param= 1.0, 0.5, 0.0, 7.0, 0.0, 0.0 => logL=   152.18517
     ## param= 1.001, 0.500, 0.000, 7.000, 0.000, 0.000 => logL=   153.21439
     ## param= 0.999, 0.500, 0.000, 7.000, 0.000, 0.000 => logL=        -Inf
     ## Error in optim(start, logL, lower = lower, upper = upper, method = optim.method,  :
     ##   non-finite finite-difference value [1]
     ## Calls: system.time ... fitCopula -> .local -> <Anonymous> -> fitCopula.ml -> optim

   ## MM: Note that 0.999 is not legal for a gumbelCopula (-> gives Error --> -Inf )
    coef(f.w); logLik(f.w)
    summary(f.w)
    stopifnot(exprs = {
        all.equal(  coef(f.),   coef(f.w), tol=0.0006) # seen 0.000187
        all.equal(logLik(f.), logLik(f.w), tol= 4e-8)  # seen 1.89e-9
    })
})

if(doExtras) withAutoprint({
    stG <- system.time(
        fGG <- fitCopula(mGG, uGG, method = "ml", estimate.variance=TRUE, traceOpt=1))
    stG # 3.7 (lynne, 2020)
    coef(fGG) ; logLik(fGG)
    stopifnot(exprs = {
        all.equal( ## dput(signif(*, 8)):
            c(m1.alpha = 2.0551924, m2.alpha = 1.5838548, w1 = 0.50994422, w2 = 0.49005578),
            coef(fGG), tol = 1e-7)
        all.equal(structure(256.799629, nobs = 1000L, df = 4L, class = "logLik"),
                  logLik(fGG), tol = 8e-8) # 9.67 e-10)
    })
    ## these two now work with "cheating" warning :
    summary(fGG)
    vcov(fGG)
    ## these now use 'l' aka 'lambda'-space :
    summary(fGG, orig=FALSE)
    (vcov(fGG, orig=FALSE) -> vGG) # "l-space"
    cov2cor(vGG)
    coef(fGG, orig=FALSE) # "l-space"
})


if(doExtras) withAutoprint({ # slowish
    st1 <- system.time(
        ff <- fitCopula(mC, uM, optim.method = "Nelder-Mead", optim.control=optCtrl))
    ## converges (vcov : Error in solve(..)... (rec.cond.number 4.783e-17)
    st1 # 11 sec (then lynne, 2020: 4.8'')
    logLik(ff)
    summary(ff)
    ## now using 'mC.' : -----
    system.time(
        ff.w <- fitCopula(mC., uM, optim.method = "Nelder-Mead", optim.control=optCtrl))
    coef(ff.w) ; logLik(ff.w)
    summary(ff.w)
    stopifnot(exprs = {
        all.equal(  coef(ff),   coef(ff.w), tol=0.008) # seen 0.00186
        all.equal(logLik(ff), logLik(ff.w), tol= 4e-8)  # seen 7.9e-9
        })
})


if(doExtras) withAutoprint({ # slowish
    system.time( # 28 sec
        f2 <- fitCopula(mC, uM, optim.method = "L-BFGS-B",    optim.control=optCtrl))
    ## converges (vcov: Error in solve(..)... (rec.cond.number 3.691e-17)
    coef(f2) ; coef(f2, orig=FALSE); logLik(f2)
    summary(f2)

    ## now using 'mC.' : -----
    system.time(
        f2.w <- fitCopula(mC., uM, optim.method = "Nelder-Mead", optim.control=optCtrl))
    coef(f2.w) ; coef(f2.w, orig=FALSE) ; logLik(f2.w)
    summary(f2.w)
    stopifnot(exprs = { ## different, f2 was "poor" :
        all.equal(  coef(f2),   coef(f2.w), tol=0.004) # seen 0.00159
        all.equal(coef(f2  , orig=FALSE),
                  coef(f2.w, orig=FALSE),   tol=0.004) # 0.001639
        all.equal(logLik(f2), logLik(f2.w), tol=0.002) # seen 2.26e-7
    })
})


## === Partially fixed parameters -- the hard cases  ===================================

RNGversion("4.0.0") # back to normal

(tX4 <- tCopula( 0.2, df = 5, df.fixed=TRUE))
(tn3 <- tCopula(-0.5, df = 3))
getTheta(tX4, attr = TRUE) # freeOnly = TRUE is default
## --> *not* showing df=5 as it is fixed
(m3 <- mixCopula(list(normalCopula(0.4), tX4, tn3), w = (1:3)/6))
                                        # -> shows 'm2.df := 5' as fixed !
(th. <- getTheta(m3, attr = TRUE))# ditto
(th  <- getTheta(m3, named= TRUE))
##' an inverse function of which(.) :
trueAt <- function(i, n) { r <- logical(n); r[i] <- TRUE; r }
## Functionality check of trueAt() :
set.seed(17); summary(Ns <- rpois(1000,3))
table(Ns) ; M <- max(Ns)
for(n in Ns) {
    i <- sort(sample.int(M, n))
    stopifnot(identical(i, which(trueAt(i, M))))
}# takes ~ 0.05 sec

stopifnot(exprs = {
    identical(th, structure(th., param.lowbnd=NULL, param.upbnd=NULL))
    identical(names(th), c("m1.rho.1", "m2.rho.1", "m3.rho.1", "m3.df",
                           paste0("w", 1:3)))
    identical(isFree(m3), !trueAt(3, n=8))# free everywhere but at [3]
})

set.seed(47)
Ut <- rCopula(2^12, m3)
##
fixedParam(m3) <- trueAt(c(3, 8), n=8)
stopifnot(exprs = {
    nParam(m3, freeOnly =FALSE) == 8
    nParam(m3, freeOnly = TRUE) == 6
    length(print( getTheta(m3, freeOnly=FALSE, named=TRUE))) == 8
    length(print( getTheta(m3, freeOnly=TRUE,  named=TRUE))) == 6
})
(iniP <- getIniParam(m3, Ut, default=NA)) # matching the freeOnly case:
stopifnot(exprs = {
    identical(iniP, getIniParam(m3, Ut))
    length(iniP) == 6
    iniP == c(rep(iniP[1], 3), 3, 1/4, 1/4)
})
## does the inital value possibly "work" ?
(llIni <- loglikCopula(u=Ut, copula=setTheta(m3, iniP, na.ok=FALSE))) # 205.0475; yes, good, ...
if(doExtras) withAutoprint({ ## try with the true parameters of 'm3'
    system.time( # 28 sec
        f38 <- fitCopula(m3, Ut, traceOpt=TRUE) )
    coef(f38) ; logLik(f38)
    summary(f38)

    ## More realistically, starting from our simple iniP[]  (instead from the true \theta !!)
    system.time(
    ffI <- fitCopula(m3, Ut, start = iniP, optim.method = "BFGS",
                     optim.control=list(maxit=999), traceOpt=TRUE) )
    coef(ffI) ; logLik(ffI)
    summary(ffI)
})


fixedParam(m3) <- trueAt(c(3, 7), n=8) # (fixed wgt in the middle)
## works with fixed at c(3,6) ...
(iniP37 <- getIniParam(m3, Ut))
if(doExtras) withAutoprint({ # *much* slower  (less parameters, why ??)
    system.time( # 28 sec
        f37 <- fitCopula(m3, Ut, start=iniP37, traceOpt=TRUE) )
    coef(f37) ; logLik(f37)
    summary(f37)
})

## setting "arbitrary" (pretty wrong) initial values :
mm <- setTheta(m3, c(0, 0, 0, df=13, w1 = .6, w2 = .4))
## (not yet ... getIniParam() is more relevant !)
