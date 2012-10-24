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

(doExtras <- copula:::doExtras())

## From source(system.file("test-tools-1.R", package = "Matrix")) :
showProc.time <- local({
    pct <- proc.time()
    function() { ## CPU elapsed __since last called__
	ot <- pct ; pct <<- proc.time()
	cat('Time elapsed: ', (pct - ot)[1:3],'\n')
    }
})


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
	 (f2.r <- fitCopula(tCopula(), uu, method="irho"))#
if(doExtras) {
    print(f2.m <- fitCopula(tCopula(), uu, method=  "ml"))# gives SE for 'df' {from optim()}
    print(f2.M <- fitCopula(tCopula(), uu, method= "mpl"))# no SE for 'df' (for now ..)
}
showProc.time()

## d = 3 : -------------
## ok with df.fixed
tC3f <- tCopula(c(.2,.7, .8), dim=3, dispstr="un", df.fixed=TRUE)
print(f3 <- fitCopula(tC3f, u3, method="itau"))

tC3 <- tCopula(c(.2,.7, .8), dim=3, dispstr="un")
	 (f3.t <- fitCopula(tC3, u3, method="itau"))
	 (f3.r <- fitCopula(tC3, u3, method="irho"))
if(doExtras) {
    print(f3.m <- fitCopula(tC3, u3, method=  "ml"))
    print(f3.M <- fitCopula(tC3, u3, method= "mpl"))
}

showProc.time()

if(!doExtras && !interactive()) q(save="no") ## so the following auto prints
##--------------------------------------------------------------------------

## d = 2 :
try( ## fails for tau = 0.8 in optim(), "non-finite finite-difference" ... FIXME
rtx <- tstFit1cop(tCopula(df.fixed=TRUE), tau.set=c(.4, .8), n.set=c(10, 25), N=64)
)
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
fitCopula(tevCopula(),		    x)
fitCopula(tevCopula(), 		    x, estimate.variance=FALSE)
try(
fitCopula(tevCopula(), 		    x, method="ml")
)

set.seed(7)
try(
rtevx <- tstFit1cop(tevCopula(, df.fixed=TRUE),
                    tau.set= c(.5, .75), n.set=c(10, 25), N=32)
)##--> singular linear system (Lapack ...)
## now "non-finite finite-difference {in optim()}

## for df.fixed=FALSE, have 2 parameters ==> cannot use "fit1":
## ....
## .... TODO

showProc.time()


set.seed(121)

### Fitting  multivariate incl margins --- mvdc
gumbelC <- gumbelCopula(3, dim=2)
gMvGam <- mvdc(gumbelC, c("gamma","gamma"), param = list(list(2,3), list(4,1)))
gMvGam # now nicely show()s -- the *AUTO*-constructed parameter names
X <- rMvdc(16000, gMvGam)
plot(X, cex = 1/4)

persp  (gMvGam, dMvdc, xlim = c(0,4), ylim=c(0,8)) ## almost discrete ????
contour(gMvGam, dMvdc, xlim = c(0,2), ylim=c(0,8))
points(X, cex = 1/16, col=adjustcolor("blue", 0.5))


if(FALSE)# unfinished
fMv <- fitMvdc(X, gMvGam)


showProc.time()
