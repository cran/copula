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
source(system.file("Rsource", "tstFit-fn.R", package="copula", mustWork=TRUE))
source(system.file("Rsource", "utils.R",     package="copula", mustWork=TRUE))
##-> assertError(), assert.EQ(), ... showProc.time()

(doExtras <- copula:::doExtras())


uu <- array(c(9, 7, 8, 3, 2,   4, 1, 5, 6, 10,
              6, 9, 1, 7, 3,   2, 5, 8, 4, 10), dim = c(10L, 2L)) / 11
set.seed(7)
u3 <- cbind(uu, round(runif(10),2))

### t-copula instead of normal -- minimal set for testing here:
## d = 2
(f1 <- fit1(tCopula(df.fixed=TRUE), x = uu))
stopifnot(identical(f1, fit1(tCopula(df.fixed=TRUE),
			     x = data.frame(uu))))# *WITH* a warning
## did not work with data.frame before 2012-08-12

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
(f2.t <- fitCopula(tCopula(), uu, method="itau"))#
if(FALSE) ## 14.Jan.2016: irho(<tCopula>) has been "non-sense" -> now erronous:
    f2.r <- fitCopula(tCopula(), uu, method="irho")
if(doExtras) {
    print(f2.m <- fitCopula(tCopula(), uu, method=  "ml"))# gives SE for 'df' {from optim()}
    print(f2.M <- fitCopula(tCopula(), uu, method= "mpl"))# no SE for 'df' (for now ..)
}
showProc.time()

## d = 3 : -------------
## ok with df.fixed
tC3f <- tCopula(c(.2,.7, .8), dim=3, dispstr="un", df.fixed=TRUE)
tC3  <- tCopula(c(.2,.7, .8), dim=3, dispstr="un")
print(f3 <- fitCopula(tC3f, u3, method="itau"))
(f3.t <- fitCopula(tC3, u3, method="itau"))# warning: coercing to df.fixed=TRUE
if(FALSE) ## irho(<tCopula>) now an error:
(f3.r  <- fitCopula(tC3, u3, method="irho"))

if(doExtras) {
    print(f3.m <- fitCopula(tC3, u3, method=  "ml"))
    print(f3.M <- fitCopula(tC3, u3, method= "mpl"))
}

showProc.time()

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
## rather need "un" dispersion: Now with smarter tCopula():
(uc4 <- tCopula(dim=d, df=nu, disp = "un", df.fixed=FALSE))
validObject(uc4p <- setTheta(uc4, value = c(P, df=nu)))
U. <- pobs(rCopula(n=1000, copula=uc4p))
splom2(U.) # => now correct dependency
(cU <- cor(U., method="kendall")) # => correct:
stopifnot(cor(P, cU[lower.tri(cU)]) > 0.99)

## Fitting a t-copula with "itau.mpl" with disp="un"
(fm4u <- fitCopula(uc4, U., method="itau.mpl", traceOpt = TRUE))
if(FALSE) { # The following are not available (yet); see ~/R/fitCopula.R
    ## Fitting a t-copula with "itau.mpl" with disp="ex"
    uc4.ex <- tCopula(dim=d, df=nu, disp = "ex", df.fixed=FALSE)
    validObject(uc4p.ex <- setTheta(uc4.ex, value = c(0.75, df=nu)))
    U.ex <- pobs(rCopula(n=1000, copula=uc4p.ex))
    (fm4e <- fitCopula(uc4.ex, U.ex, method="itau.mpl"))
    ## Fitting a t-copula with "itau.mpl" with disp="ar"
    uc4.ar <- tCopula(dim=d, df=nu, disp = "ar1", df.fixed=FALSE)
    validObject(uc4p.ar <- setTheta(uc4.ar, value = c(0.75, df=nu)))
    U.ar <- pobs(rCopula(n=1000, copula=uc4p.ar))
    (fm4e <- fitCopula(uc4.ar, U.ar, method="itau.mpl"))
}

## Extra checks --------------------------------------------------------
if(doExtras) { ## Typically not run on R CMD check
tCop <- tCopula(c(0.2,0.4,0.6), dim=3, dispstr="un", df=5)
set.seed(101)
x <- rCopula(n=200, tCop) # "true" observations (simulated)
## Maximum likelihood (start = (rho[1:3], df))
print(tc.ml  <- fitCopula(tCopula(dim=3, dispstr="un"), x, method="ml",
                          start=c(0,0,0, 10)))
print(tc.ml. <- fitCopula(tCopula(dim=3, dispstr="un"),
                          x, method="ml")) # without 'start'
## Maximum pseudo-likelihood (the asymptotic variance cannot be estimated)
u <- pobs(x)
print(tc.mpl <- fitCopula(tCopula(dim=3, dispstr="un"),
                          u, method="mpl", estimate.variance=FALSE,
                          start= c(0,0,0,10)))
## Without 'start'
tc.mp. <- fitCopula(tCopula(dim=3, dispstr="un"), x, estimate.variance=FALSE)
print(tc.mp.)
noC <- function(x) { x@fitting.stats$counts <- NULL ; x }

assert.EQ(noC(tc.ml) , noC(tc.ml.), tol= .005)
assert.EQ(noC(tc.mpl), noC(tc.mp.), tol= .100, giveRE=TRUE)

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
} # end Xtras




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

if(FALSE)# unfinished --- TODO maybe move below ('doExtras')!
fMv <- fitMvdc(X, gMvGam) # needs 'start' ...

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
)
## must not valid: stopifnot(!isTRUE(validObject(mvW, test=TRUE)))

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
par <- 2^c((0:32)/16, 2+(1:10)/8)
llg <- sapply(par, function(p) loglikCopula(param=p, u=u, copula=gumbel))
if(dev.interactive()) plot(par, llg, type="b", col=2)
stopifnot(diff(llg) < 0) # so the maximum is for par = 2^0 = 1 --> at *boundary* of interval
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
if(FALSE) ## takes ??? (> one hour)
## try  maxiter =  / trace =  to see what's happening
system.time(fit <- fitCopula(ellipCopula("t", dim=d, dispstr="un"),
                             data=U, method="mpl"))


## > Warning in makePosDef(sigma, delta = 0.001) :
## >   Estimate is not positive-definite. Correction applied.
## > Process R killed: 9 at Sat Nov 14 19:48:55 2015
1
## ... so R crashes (this also happens if the data is much 'closer' to a
## t-copula than my data above, so that's not the problem). I'm wondering
## about the following.

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
(tevc <- tevCopula(iTau(tevCopula(), 0.75)))
set.seed(1); str(x <- rCopula(1000, tevc))
plot(x, main = "1000 samples of tevCopula(iTau(tevCopula(), 0.75))")
fitCopula(tevCopula(),		    x, method="irho")# warning
fitCopula(tevCopula(df.fixed=TRUE), x, method="itau")# fine
fitCopula(tevCopula(df.fixed=TRUE), x)# two warnings ==> do not estimate.var:
fitCopula(tevCopula(df.fixed=TRUE), x, estimate.variance=FALSE)
fitCopula(tevCopula(df.fixed=TRUE), x, method="ml")
## 'df' is not estimated, but it should
fitCopula(tevCopula(),		    x)# df as parameter, but "kept at df := 4" (FIXME?)
fitCopula(tevCopula(), 		    x, estimate.variance=FALSE)
fitCopula(tevCopula(), 		    x, method="ml")# currently ~= the one above

set.seed(7)
try(# now works .. slowly!
rtevx <- tstFit1cop(tevCopula(, df.fixed=TRUE),
                    tau.set= c(.5, .75), n.set=c(10, 25), N=32)
)##--> singular linear system (Lapack ...)
## now "non-finite finite-difference {in optim()}

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO
showProc.time()

data(rdj)
rdj <- rdj[,2:4]
dim(u <- pobs(rdj))# 1262 3
fc <- frankCopula(dim=3)
ffc <- fitCopula(fc, u) ## (failed in 0.999-4 {param constraints})
ffc
## summary(ffc) -- *not* useful -- FIXME: needs print.summary.fitCopula()
stopifnot(all.equal(unname(coef(ffc)),
                    2.866564929, tolerance = 1e-5))

showProc.time()

