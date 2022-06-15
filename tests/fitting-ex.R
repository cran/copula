## Copyright (C) 2021 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
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
source(system.file("Rsource", "tstFit-fn.R", package="copula", mustWork=TRUE))
source(system.file("Rsource", "utils.R",     package="copula", mustWork=TRUE))
##-> assertError(), assert.EQ(), ... showProc.time()

(doExtras <- copula:::doExtras())


uu <- array(c(9, 7, 8, 3, 2,   4, 1, 5, 6, 10,
              6, 9, 1, 7, 3,   2, 5, 8, 4, 10), dim = c(10L, 2L)) / 11
suppressWarnings(RNGversion("3.5.0")) ## << FIXME: eliminate - adapting where needed below !
set.seed(7)
u3 <- cbind(uu, round(runif(10),2))

### t-copula instead of normal -- minimal set for testing here:
## d = 2
## fit1() here catches the error: "irho" not available for "tCopula":
(f1 <- fit1(tCopula(df.fixed=TRUE), x = uu))
stopifnot(identical(f1, fit1(tCopula(df.fixed=TRUE),
			     x = data.frame(uu))))# *WITH* a warning; keep data.frame() working!
## did not work with data.frame before 2012-08-12

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
(f2.t <- fitCopula(tCopula(), uu, method="itau"))# (+ warning: ... 'df.fixed=TRUE' )
assertError( ## 14.Jan.2016: irho(<tCopula>) has been "non-sense" -> now erronous:
    fitCopula(tCopula(), uu, method="irho"), verbose=TRUE)
showProc.time()
if(doExtras) {
    print(f2.m <- fitCopula(tCopula(), uu, method=  "ml"))
    print(f2.M <- fitCopula(tCopula(), uu, method= "mpl"))
    print(summary(f2.m)) # gives SE for 'df' {from optim()}
    print(summary(f2.M)) # no SE for 'df' (for now ..)
    stopifnot(all.equal(coef(f2.m), coef(f2.M)))
    showProc.time()
}

## d = 3 : -------------
## ok with df.fixed
tC3f <- tCopula(dim=3, dispstr="un", df.fixed=TRUE)
tC3  <- tCopula(dim=3, dispstr="un")
f3   <- fitCopula(tC3f, u3, method="itau")
f3.t <- fitCopula(tC3 , u3, method="itau") # warning: coercing to df.fixed=TRUE
summary(f3.t)
cf3 <- coef(f3, SE = TRUE)
stopifnot(all.equal(unname(cf3), cbind(c(0.374607, 0.309017, 0.374607),
                                       c(0.386705, 0.325995, 0.405493)),
                    tol = 5e-5), # seen 6e-7
          all.equal(coef(f3), coef(f3.t)))
if(FALSE) ## Error: iRho() method for class "tCopula" not yet implemented
(f3.r  <- fitCopula(tC3, u3, method="irho"))
showProc.time()

octr <- list(trace=TRUE, REPORT=1, reltol = 1e-5)# less accurate than by default
if(doExtras) {
    fitC.m <- function(C, data, method, ...)
        fitCopula(C, data, method=method, optim.control=octr, traceOpt=TRUE, ...)
    print(f3.m <- fitC.m(tC3, u3, method=  "ml")); c.m <- coef(f3.m, SE=TRUE)
    print(f3.M <- fitC.m(tC3, u3, method= "mpl")); c.M <- coef(f3.M, SE=TRUE)
    showProc.time()
    is.matrix(print(c.m))
    is.character(SE <- "Std. Error")
    stopifnot(all.equal(c.m[,1], c.M[,1]), # the estimates don't differ; the SE's do:
              c.M[-4,SE] >= c.m[-4,SE],
              c.m[, SE]  > 0)
}

set.seed(17)
d <- 5 # dimension
nu <- 4 # degrees of freedom
## define and sample the copula, build pseudo-observations
ec4 <- tCopula(dim=d, df=nu, df.fixed=TRUE) # <- copula with param NA
(r <- iTau(ec4, tau <- c(0.2, 0.4, 0.6)))
P <- c(r[2], r[1], r[1], r[1], # upper triangle (w/o diagonal) of corr.matrix
             r[1], r[1], r[1],
                   r[3], r[3],
                         r[3])
assertError( setTheta(ec4, value = P) )
validObject(ex4. <- setTheta(ec4, value = 0.8)) # set the only (non-fixed) parameter
## TODO "check" getTheta(ex4., ...)
## rather need "un" dispersion: Now with smarter tCopula():
(uc4 <- tCopula(dim=d, df=nu, disp = "un", df.fixed=FALSE))
validObject(uc4p <- setTheta(uc4, value = c(P, df=nu)))
U. <- pobs(rCopula(n=1000, copula=uc4p))
splom2(U.) # => now correct dependency
(cU <- cor(U., method="kendall")) # => correct:
stopifnot(cor(P, cU[lower.tri(cU)]) > 0.99)

## Fitting a t-copula with "itau.mpl" with disp="un"
(fm4u <- fitCopula(uc4, U., method="itau.mpl", traceOpt = TRUE))
## Fitting  t-copulas  .............  with disp = "ex" and "ar" :
uc4.ex <- tCopula(dim=d, df=nu, disp = "ex", df.fixed=FALSE)
uc4.ar <- tCopula(dim=d, df=nu, disp = "ar1", df.fixed=FALSE)
validObject(uc4p.ex <- setTheta(uc4.ex, value = c(0.75, df=nu)))
validObject(uc4p.ar <- setTheta(uc4.ar, value = c(0.75, df=nu)))
U.ex <- pobs(rCopula(n=1000, copula=uc4p.ex))
U.ar <- pobs(rCopula(n=1000, copula=uc4p.ar))
if(FALSE) { # The following are not available (yet); see ~/R/fitCopula.R
    ## Fitting a t-copula with "itau.mpl" with disp="ex"
    (fm4e <- fitCopula(uc4.ex, U.ex, method="itau.mpl"))
    ## Fitting a t-copula with "itau.mpl" with disp="ar"
    (fm4e <- fitCopula(uc4.ar, U.ar, method="itau.mpl"))
}

## Extra checks --------------------------------------------------------
if(doExtras) { ## Typically not run on R CMD check
tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5)
set.seed(101)
x <- rCopula(n=200, tCop) # "true" observations (simulated)
## Maximum likelihood (start = (rho[1:3], df))
print(summary(tc.ml <-
                  fitCopula(tCopula(dim=3, dispstr="un"), x, method="ml",
                            start=c(0,0,0, 10))))
print(summary(tc.ml. <-
                  fitCopula(tCopula(dim=3, dispstr="un"), x, method="ml")))# w/o 'start'

## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
u <- pobs(x)
print(tc.mpl <- fitCopula(tCopula(dim=3, dispstr="un"),
                          u, method="mpl", estimate.variance=FALSE,
                          start= c(0,0,0,10)))
## Without 'start'
tc.mp. <- fitCopula(tCopula(dim=3, dispstr="un"), x, estimate.variance=FALSE)
print(tc.mp.)
noC <- function(x) { x@fitting.stats$counts <- NULL; x@call <- quote(dummy()); x }

assert.EQ(noC(tc.ml) , noC(tc.ml.), tol= .005) # nothing
assert.EQ(noC(tc.mpl), noC(tc.mp.), tol= .100, giveRE=TRUE) # shows diff

## The same t copula but with df.fixed=TRUE (=> use the same data!)
tC3u5 <- tCopula(dim=3, dispstr="un", df=5, df.fixed=TRUE)
## Maximum likelihood (start = rho[1:3])
print(tcF.ml  <- fitCopula(tC3u5, x, method="ml", start=c(0,0,0)))
print(tcF.ml. <- fitCopula(tC3u5, x, method="ml"))  # without 'start'
assert.EQ(noC(tcF.ml), noC(tcF.ml.), tol= 4e-4)
print(vcov(tcF.ml)) # the (estimated, asymptotic) var-cov matrix

## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
print(tcF.mpl <- fitCopula(tC3u5, u, method="mpl", estimate.variance=FALSE,
                           start=c(0,0,0)))
print(tcF.mp. <- fitCopula(tC3u5, u, method="mpl", estimate.variance=FALSE))
assert.EQ(noC(tcF.mpl), noC(tcF.mp.), tol = 1e-5)
} # end Extras ----------------------------------------------------

## fitMvdc() -- first 2 D -- from Yiyun Shou's bug report: ---------------------

ct.2 <- tCopula(param=0.2, dim=2, dispstr = "ex", df.fixed=TRUE)
mvt.2.ne <- mvdc(copula = ct.2, margins = c("norm", "exp"),
                 paramMargins = list(list(mean = 0, sd = 2), list(rate = 2)))
mvt.2.ne ## --> four free parameters in total: rho, mean, sd, and rate:

getTheta(mvt.2.ne, attr = TRUE) # - shows only the copula parameter -- FIXME !

## simulate data and fit:
set.seed(17); x.samp <- rMvdc(250, mvt.2.ne)
fit2ne <- fitMvdc(x.samp, mvt.2.ne, start= c(1,1,1, rho = 0.5),
                  optim.control = list(trace = TRUE), hideWarnings=FALSE)
summary(fit2ne)
(confint(fit2ne) -> ci.2ne)
stopifnot(exprs = {
    all.equal(coef(fit2ne),
              c(m1.mean=0.061359521, m1.sd=2.0769423,
                m2.rate=2.0437937, rho.1 = 0.15074002), tol=1e-7)# seen 1.48e-8
    all.equal(c(ci.2ne),
              c(-0.18309, 1.9064, 1.8019, 0.012286,
                 0.30581, 2.2474, 2.2857, 0.28919), tol = 4e-4) # seen 1.65e-5
})



### Fitting  multivariate incl margins --- mvdc --------------------------------------
### ===========================================

set.seed(121)
gumbelC <- gumbelCopula(3, dim=2)
gMvGam <- mvdc(gumbelC, c("gamma","gamma"), param = list(list(2,3), list(4,1)))
gMvGam # now nicely show()s -- the *AUTO*-constructed parameter names
stopifnot(identical(gMvGam@paramMargins,
                    list(list(shape = 2, rate = 3),
                         list(shape = 4, rate = 1))))
X <- rMvdc(16000, gMvGam)
smoothScatter(X, main = "rMvdc(1600, gMvGam)")

persp  (gMvGam, dMvdc, xlim = c(0,4), ylim=c(0,8)) ## almost discrete ????
contour(gMvGam, dMvdc, xlim = c(0,2), ylim=c(0,8))
points(X, pch = ".", cex = 2, col=adjustcolor("blue", 0.5))

if(doExtras) {
    st <- system.time(
        fMv <- fitMvdc(X, gMvGam, start = c(1,1,1,1, 1.3),# method="BFGS",
                       optim.control= list(trace=TRUE)))
    print(st) # ~ 59 sec. (lynne 2015); ~ 47 sec (2020)
    print(summary(fMv))
    cfM <- coef(fMv, SE=TRUE)
    print(    all.equal(c(2,3, 4,1, 3), unname(cfM[,"Estimate"]))) # 0.01518
    stopifnot(all.equal(c(2,3, 4,1, 3), unname(cfM[,"Estimate"]), tol = 0.02))
}

pFoo <- function(x, lower.tail=TRUE, log.p=FALSE)
     pnorm((x - 5)/20, lower.tail=lower.tail, log.p=log.p)
dFoo <- function(x, lower.tail=TRUE, log.p=FALSE)
     1/20* dnorm((x - 5)/20, lower.tail=lower.tail, log.p=log.p)
qFoo <- qunif # must exist; not used for fitting

## 'Foo' distribution has *no* parameters:
mv1 <- mvdc(gumbelC, c("gamma","Foo"), param= list(list(3,1), list()))
validObject(mv1)
stopifnot(nrow(R <- rMvdc(3, mv1)) == 3, ncol(R) == 2)
## a wrong way:
assertError(
  mvW <- mvdc(gumbelC, c("gamma","Foo"), param= list(list(3,1), list(NULL)))
, verbose=TRUE)
## >> invalid class “mvdc” object: 'paramMargins[[2]]' must be named properly


showProc.time()

## An example which fails (and should)? --
## From: Suzanne Li  @...queensu.ca, Date: Fri, 9 Aug 2013
gumbel <- archmCopula(family = "gumbel",dim = 2)

set.seed(47)# MM {not important, but still want sure reproducibility}
u <- cbind(runif(10),runif(10))
cor(u[,1], u[,2], method="kendall")
## [1] -0.02222222 -- slightly negative
## this now gives an error:
try(fgu <- fitCopula(gumbel, u, method = "ml"))
## Error in optim(start, loglikCopula, lower = lower, upper = upper, method = method,  :
##   non-finite finite-difference value [1]
## In addition: Warning message:
## In .local(copula, tau, ...) : tau is out of the range [0, 1]
copGumbel@paraInterval # -> [1, Inf) = exp([0, Inf))
length(par <- 2^c((0:32)/16, 2+(1:10)/8)) # 43
(t1 <- system.time(
llg <- sapply(par, function(p) loglikCopula(param=p, u=u, copula=gumbel))
))
(t2 <- system.time(
llg. <- loglikCopulaMany(par, u=u, copula=gumbel)
))
all.equal(llg, llg., tol=0)
if(dev.interactive()) plot(par, llg, type="b", col=2)
stopifnot(exprs = {
    diff(llg) < 0 # so the maximum is for par = 2^0 = 1 --> at *boundary* of interval
    all.equal(llg, llg., tol=4e-16) # see 0 diff (Lnx 64b)
})

## FIXME -- "ml" should return the boundary case, or a much better error message
## These work (with a warning {which is interesting, but maybe should not be a warning}
## "perfectly": They give the correct boundary case:
fg.itau <- fitCopula(gumbel, u, method = "itau")
fg.irho <- fitCopula(gumbel, u, method = "irho")

## Is it just this problem?
## well, the likelihood was not ok for large param=theta; now is:
lrgP <- 100*seq(8, 100, by = 3)
llrg <- vapply(lrgP, function(P) loglikCopula(param=P, u=u, copula=gumbel), NA_real_)
stopifnot(is.finite(llrg), diff(llrg) < 0, llrg < -11990)## no longer NaN
if(FALSE)
    plot(lrgP, llrg, type="b", col=2) ## all fine

## Now is it because we *really* should use  elme()  and the "nacopula" families?
## No, this fails too: "outside interval" error {instead of just -Inf !}:
(copG <- onacopulaL("Gumbel", list(NaN,1:2)))
## Estimation -> error for now -- (FIXME!)
try(efm <- emle(u, copG))


## A simple example with *negative* correlation;
## Want *more helpful* error message here:
set.seed(7)
u1 <- seq(0,1, by=1/128)[2:127]; u2 <- -u1 + round(rnorm(u1)/4,2); u <- pobs(cbind(u1,u2))
plot(u)
msg <- tryCatch(fitCopula(gumbelCopula(), data = u), error=function(e)e$message)
## check error message __FIXME__ want "negative correlation not possible"
## or NO ERROR and a best fit to tau=0 [and the same for other Archimedean families!]
msg
showProc.time()


## Date: Sat, 14 Nov 2015 20:21:27 -0500
## Subject: Replace makePosDef() by nearPD()?
## From: Marius Hofert <marius.hofert@uwaterloo.ca>
## To: Ivan Kojadinovic <ivan.kojadinovic@univ-pau.fr>, Jun Yan
## 	<jun.yan@uconn.edu>, Martin Maechler <maechler@stat.math.ethz.ch>

## I currently face the problem of estimating a t-copula to a (504, 413)
## data matrix. A MWE is:

require(copula)
set.seed(271)
n <- 504
d <- 413
d <- 50 # more realistic -- should become faster !!
U <- matrix(runif(n*d), nrow=n, ncol=d)
if(FALSE) ## takes ??? (> one hour); *control = octr : some trace, not too accurate
system.time(fit <- fitCopula(ellipCopula("t", dim=d, dispstr="un"),
                             data=U, method="mpl", optim.control = octr))

## > Warning in makePosDef(sigma, delta = 0.001) :
## >   Estimate is not positive-definite. Correction applied.
## > Process R killed: 9 at Sat Nov 14 19:48:55 2015
1
## ... so R crashes (this also happens if the data is much 'closer' to a
## t-copula than my data above, so that's not the problem). I'm wondering
## about the following.
showProc.time()

## 4) Implement: first estimating (via pairwise inversion of
## Kendall's tau) the dispersion matrix and then estimating the d.o.f.
## via MLE (based on the fitted dispersion matrix). Already in use in
## vignettes --- would be good to have this in the
## package as. What we could do is to extend
## fitCopula.icor(, method="kendall"): if the copula has a d.o.f. parameter, then don't
## treat it as fixed but estimated it via MLE... (or to get this behavior
## via a method argument or so?)

## 6) After 4), lambda() returns a vector of length 2* d(d-1)/2 ... ok
##    in the bivariate case but not in higher dimensions => want a list

###----------- xvCopula()  [ <--> ../man/xvCopula.Rd ] ---------
set.seed(12)
x <- rCopula(64, claytonCopula(pi))# "small" n ==> fast
fG <- fitCopula(gumbelCopula(), x)
(v <- c(logL   = logLik(fG),
	xv.loo = xvCopula(gumbelCopula(), x       ), # l.o.o. CV
	xv.5   = xvCopula(gumbelCopula(), x, k = 5)))# 5-fold CV
stopifnot(all.equal(unname(v), c(32.783677, 32.835744, 32.247463),
		    tolerance = 1e-7))

## now also with rotCopula: --
## NB: upper=NA / upper=Inf is *better* for optim(*, "L-BFGS-B")  than upper= "Large" !!!
## --- (this is a puzzle in optim() itself !)
if(FALSE)# has been useful:
trace(copula:::fitCopula.ml,
      trace = quote({cat("in fitCopula.ml(): copula=\n");print(copula)})
      )
if(FALSE) ## explicit upper=Inf "works" (*not* finding very good value:
    print(ff <- fitCopula(rotCopula(claytonCopula()), upper=NA, x, traceOpt=TRUE))
if(FALSE) ## explicit upper=Inf or =NA work; upper=100 "fails" giving bad result !
    print(f1c <- fitCopula(rotCopula(claytonCopula()), upper=100, x, traceOpt=TRUE))
if(FALSE) ## "Brent" instead of default "L-BFGS-B" -- "works" very nicely with upper=100,
    print(fB <- fitCopula(rotCopula(claytonCopula()), upper=100, x,
                          optim.method="Brent",traceOpt=TRUE))
## but *not* with upper=Inf or NA (--> error: .. must be finite..) or upper = 1e300  !!
try( fitCopula(rotCopula(claytonCopula()), upper=Inf, x, optim.method="Brent"))
try( fitCopula(rotCopula(claytonCopula()), upper=NA,  x, optim.method="Brent"))
if(FALSE) ## very bad ('start' is not used with "Brent", but the upper bound is!)
    print(fBL <- fitCopula(rotCopula(claytonCopula()), upper=1e300, x,
                           optim.method="Brent",traceOpt=TRUE))

## this "works" with 'Inf', but not with default --> upper = <large> !
(xvR <- xvCopula(rotCopula(claytonCopula()), upper = Inf, x, k = 8, verbose=FALSE))
stopifnot(all.equal(xvR, 21.050569, tolerance = 1e-6))

## Now 'copula2 = indepCopula(..) gets 'dim = 3' from the first:
kh.C <- khoudrajiCopula(claytonCopula(2.5, dim=3), shapes = (1:3)/4)
kh.C
getTheta(kh.C, attr=TRUE) # bounds look good
x3 <- cbind(x, x[,1]^2 + runif(nrow(x)))
u3 <- pobs(x3)
if(FALSE) # Error:  "non-finite finite-difference [3]"
xvK <- xvCopula(kh.C, u3, k = 8)
if(FALSE) # fails (same Error)
fK <- fitCopula(kh.C, u3)
## "works" :
p1 <- (1:20)/2
l1 <- sapply(p1, function(th1) loglikCopula(c(th1, (1:3)/4), u3, kh.C))
plot(p1, l1, type = "b") # nice maximum at around 6.xx
showProc.time()

## works two:
fixedParam(kh.C) <- c(FALSE, FALSE, TRUE,TRUE)
getTheta(kh.C)
summary(fK12 <- fitCopula(kh.C, u3))
fixedParam(kh.C) <- FALSE # all free now
fK4 <- fitCopula(kh.C, u3, start = c(coef(fK12), 0.3, 0.5),
                 optim.method = "L-BFGS-B", traceOpt=TRUE)
summary(fK4)
## -> shape1 ~= shape2 ~= 0
## dput(signif(unname(coef(fK4)), 8))
c4 <- c(4.4482535, 0, 0, 0.63730128)
all.equal(c4, unname(coef(fK4)), tol=0)# 7.342e-9
stopifnot(all.equal(c4, unname(coef(fK4)), tol=4e-8))
showProc.time()

kh.r.C <- khoudrajiCopula(rotCopula(claytonCopula(1.5)), shapes = c(1,3)/4)
u.krC <- rCopula(999, kh.r.C)
fm.krC <- fitCopula(copula = kh.r.C, data = u.krC, traceOpt=TRUE)# works (w/ warning)
xvK.R  <- xvCopula(kh.r.C, u.krC, k = 10)
showProc.time()


## From: Ivan, 27 Jul 2016 08:58
set.seed(123)
u <- pobs(rCopula(300, joeCopula(4)))
fjc <- fitCopula(joeCopula(), data = u, method = "itau")
## Now a warning.  Previously gave 'Error in dCor(cop) :'
##   dTau() method for class "joeCopula" not yet implemented
summary(fjc)
all.equal(target = c(alpha = 4),  coef(fjc))
fjcM <- fitCopula(joeCopula(), data = u, start = coef(fjc))
confint(fjcM)# 2.5 %   97.5 %
##    alpha 3.077266 4.294797
summary(fjcM)
if(FALSE)
 dput(signif(drop(coef(fjcM, SE=TRUE)), 4))
stopifnot(
    all.equal(c(Estimate = 3.686, `Std. Error` = 0.31),
              drop(coef(fjcM, SE=TRUE)), tol = 0.0002) # seen 0.000158
)


if(!doExtras) q(save="no") ## so the following auto prints
##--------------------------------------------------------------------------
## TODO: no longer use q() / quit() -- because it breaks 'covr'
## --- R-devel Q: How can we turn on auto-printing inside an if() { .. } clause via  withVisible()?


## d = 2 :
## ----- catching fitCopula() error: 'Lapack routine dgesv: system is exactly singular: U[2,2] = 0'
##_FIXME_ since ca.end of 2015, this gives *different error*: iRho() not available for 'tCopula's'
rtx <- tstFit1cop(tCopula(df.fixed=TRUE), tau.set=c(.4, .8),
                  n.set= c(10, 25), N=64)
## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO

showProc.time()

## The other example  with 'df' (fixed / free):
(tevc <- tevCopula(iTau(tevCopula(), 0.75))) ##   rho.1 = 0.9751784, df = 4
set.seed(1); str(x <- rCopula(1000, tevc))
plot(x, main = "1000 samples of tevCopula(iTau(tevCopula(), 0.75))")
mirho <- fitCopula(tevCopula(),		     x, method="irho")# warning -> 'df.fixed=TRUE"
mitau <- fitCopula(tevCopula(df.fixed=TRUE), x, method="itau")# fine
cfR <- coef(mirho, SE=TRUE)
cfT <- coef(mitau, SE=TRUE)
stopifnot(
    all.equal(0.97, cfR[[1]], tol=1e-4)# see 3.65e-5
)
system.time(
fmM. <- fitCopula(tevCopula(df.fixed=TRUE), x)# two warnings ==> do not estimate.var:
) # 23.5 sec -- the expensive part was trying get and invert the Hessian
system.time(
fmM <- fitCopula(tevCopula(df.fixed=TRUE), x, estimate.variance=FALSE)
) # 0.22 sec  !!
stopifnot( all.equal(coef(fmM), coef(fmM.), tol=1e-15)) # even tol=0
system.time(
fML <- fitCopula(tevCopula(df.fixed=TRUE), x, method="ml")
)# 0.24 sec: very fast, too
showProc.time()
## 'df' is not estimated, but it should
## see above:  tevCopula() is treated as tevCopula(df.fixed=TRUE))
system.time(
fM2  <- fitCopula(tevCopula(),		    x)# df as parameter, -- for var() kept at df = 2.53
) # 24s
system.time(
fM2. <- fitCopula(tevCopula(), 		    x, estimate.variance=FALSE)
) # 0.44 sec
system.time(
fMm <- fitCopula(tevCopula(), 		    x, method="ml")
) # 0.56
showProc.time()
c2  <- coef(fM2 , SE=TRUE)
c2. <- coef(fM2., SE=TRUE)
cm  <- coef(fMm,  SE=TRUE) ; cm # full
all.equal(c2[,1], c2.[,1], tol=0)
all.equal(c2[,1], cm [,1], tol=0)
stopifnot(exprs = {
    ## coef : 'rho'1' are the same
    all.equal(c2[,1], c2.[,1], tol=1e-15)
    all.equal(c2[,1], cm [,1], tol=1e-15)
    all.equal(c(rho.1 = 0.96206, df = 2.5306), cm[,"Estimate"], tol=5e-5) # seen 5.62e-6
    ## Std.Error do differ: dput(signif(cm[,2], 4))
    all.equal(c(rho.1 = 0.00663, df = 0.457), cm[,"Std. Error"], tol=0.002)# seen 0.0007
})

showProc.time()
set.seed(7)
system.time(# now works .. slowly!
rtevx <- tstFit1cop(tevCopula(, df.fixed=TRUE),
                    tau.set= c(.5, .75), n.set=c(10, 25), N=32)
)
## gives several "non-finite finite-difference {in optim()}
rtevx ## some NaN, NA

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO
showProc.time()

data(rdj)
rdj <- rdj[,2:4]
splom2(rdj, cex=0.4, ## this gave an error for a day or so:
       col.mat = matrix(adjustcolor("black", 0.5), nrow(rdj), 3))
dim(u <- pobs(rdj))# 1262 3
fc <- frankCopula(dim=3)
ffc <- fitCopula(fc, u) ## (failed in 0.999-4 {param constraints})
summary(ffc)
stopifnot(all.equal(2.866565, unname(coef(ffc)), tolerance = 1e-5))

showProc.time()

