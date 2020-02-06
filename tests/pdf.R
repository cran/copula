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

### pdf = dCopula()  ---  but then also pCopula()
### ===   ~~~~~~~~~                     ~~~~~~~~~

require(copula)
## library(fCopulae)

source(system.file("Rsource", "utils.R", package="copula", mustWork=TRUE))
## tryCatch.W.E() etc
## All non-virtual copula classes: ../inst/Rsource/cops.R
source(system.file("Rsource", "cops.R", package="copula", mustWork=TRUE))
## --> copcl, copObs, copBnds,  excl.2 , copO.2, copBnd.2
(doExtras <- copula:::doExtras())


### preparation for a grid
n1 <- 17
n2 <- 21
eps <- .001 ## <- not going to very extremes
u1 <- seq(0, 1, length=n1); u1[1] <- eps; u1[n1] <- 1-eps
u2 <- seq(0, 1, length=n2); u2[1] <- eps; u2[n2] <- 1-eps

### d=2 ########################################################################

Exp.grid <- function(...)
    as.matrix(expand.grid(..., KEEP.OUT.ATTRS = FALSE))

umat <- Exp.grid(u1=u1, u2=u2)

## all copulas give the same tau except amh (TODO: *range* of tau's ??)
tau <- 0.5

if(!dev.interactive(orNone=TRUE)) pdf("densCop_2d.pdf")

fCols <- colorRampPalette(c("red", "white", "blue"), space = "Lab")

options(width = 137)# -> nicer table printing
showProc.time()


## frankCopula
theta.fr <- iTau(frankCopula(), tau)
dcop <- matrix(dCopula(umat, frankCopula(param=theta.fr, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( frankCopula(%.4g) )", theta.fr))
round(dcop, 3)


## claytonCopula
(theta.cl <- iTau(claytonCopula(), tau))
stopifnot(all.equal(theta.cl, copClayton@iTau(tau), tolerance = 1e-13))
dcop <- matrix(dCopula(umat, claytonCopula(param=theta.cl, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf(" Density( claytonCopula(%.4g) )", theta.cl))
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( claytonCopula(%.4g) )", theta.cl))
round(dcop, 3)



## gumbelCopula
theta.gu <- iTau(gumbelCopula(), tau)
stopifnot(all.equal(theta.gu, copGumbel@iTau(tau), tolerance = 1e-13))
dcop <- matrix(dCopula(umat, gumbelCopula(param=theta.gu, dim = 2)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( gumbelCopula(%.4g) )", theta.gu))
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( gumbelCopula(%.4g) )", theta.gu))
round(dcop, 3)


## normalCopula
uB <- cbind(c(0:1, .5), (1:9)/10) # values at boundaries
fC <- dCopula(uB, normalCopula(0.55))
stopifnot(is.finite(fC), length(fC)==nrow(uB), fC[-3*(1:3)] == 0)
theta.n <- iTau(normalCopula(), tau)
dcop <- matrix(dCopula(umat, normalCopula(param=theta.n, dim = 2)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( normalCopula(%.4g) )", theta.n))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( normalCopula(%.4g) )", theta.n))
round(dcop, 3)


## tCopula
fC <- dCopula(uB, tCopula(0.55))
stopifnot(is.finite(fC), length(fC)==nrow(uB), fC[-3*(1:3)] == 0)
(theta.t. <- iTau(tCopula(df=10), tau))
dcop <- matrix(dCopula(umat, tCopula(param=theta.t., dim = 2, df=10)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( tCopula(%.4g, df = 10) )", theta.t.))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( tCopula(%.4g, df = 10) )", theta.t.))
round(dcop, 3)
## tCopula -- df=4
(theta <- iTau(tCopula(df=4), tau))
dcop <- matrix(dCopula(umat, tCopula(param=theta, dim = 2, df=4)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( tCopula(%.4g, df = 4) )", theta))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( tCopula(%.4g, df = 4) )", theta))
round(dcop, 3)

## galambosCopula
#
(theta <- iTau(galambosCopula(), tau))
stopifnot(all.equal(tau, tau(galambosCopula(theta)), tolerance = 1e-5))
dcop <- matrix(dCopula(umat, galambosCopula(param=theta)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( galambosCopula(%.4g) )", theta))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( galambosCopula(%.4g) )", theta))
round(dcop, 3)

## plackettCopula
#
(theta <- iTau(plackettCopula(), tau))
stopifnot(all.equal(tau, tau(plackettCopula(theta)), tolerance = 1e-5))
dcop <- matrix(dCopula(umat, plackettCopula(param=theta)),
               n1, n2)
filled.contour(u1, u2, log(dcop), color.palette = fCols,
               main=sprintf("log Density( plackettCopula(%.4g) )", theta))
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( plackettCopula(%.4g) )", theta))
round(dcop, 3)


## amhCopula
tau <- 0.3 ## -- to be in range for AMH
(theta <- iTau(amhCopula(), tau))
stopifnot(all.equal(tau, tau(amhCopula(theta)), tolerance = 1e-5))
dcop <- matrix(dCopula(umat, amhCopula(param=theta)),
               n1, n2)
filled.contour(u1, u2, dcop, color.palette = fCols,
               main=sprintf("Density( amhCopula(%.4g) )", theta))
round(dcop, 3)

showProc.time()
if(!dev.interactive()) dev.off()


### d > 2 ######################################################################

### for package fCopulae
## dcop.f <- darchmCopula(umat, alpha=2, type=1, alternative=T)
## round(matrix(dcop.f, n1, n2), 3)
## dcop.fCopulae <- darchmCopula(umat, alpha=theta.fr, type = 5, alternative = TRUE)
## round (matrix(dcop.fCopulae, n1, n2), 3)

## dim = 3, frankCopula
um3 <- Exp.grid(u1=u1, u2=u2, u3=u1)
dcE <- tryCatch(dcopula(um3), error=identity)
dcop <- dCopula(um3, frankCopula(param=theta.fr, dim = 3))
stopifnot(dcop >= 0, # no NA here
          inherits(dcE, "error"),
          grepl("defunct", conditionMessage(dcE), ignore.case=TRUE))
round(array(dcop, c(n1,n2,n1)), 3)

## dim = 4 --- fine as long we are "out of corners"
um4 <- Exp.grid(u1=u1, u2=u2, u3=u1, u4=rev(u2))
dcop <- dCopula(um4, frankCopula(param=theta.fr, dim = 4))
stopifnot(dcop >= 0,# no NA here
	  all.equal(dcop, copFrank@dacopula(um4, theta = theta.fr)))
## round(array(dcop, c(n1,n2,n1,n2)), 3)

showProc.time()


### --- now look at dCopula() and pCopula() for *all* copulas:

##' u "mostly" in [0,1] .. with exceptions that are even NA, NaN
mku <- function(n, fr = 1/20, sd = 1/10) {
    r <- runif(n) + rnorm(n, sd=sd)
    r[sample(n, ceiling(fr*n))] <- c(NA,NaN)
    r
}

## This is from ./moments.R --- keep in sync! ---
tau.s <- c(       -.1, 0, 0.05805, (1:2)/9, 0.3)
names(tau.s) <- paste0("tau=", sub("0[.]", ".", formatC(tau.s)))
suppressWarnings(
    tTau <- sapply(tau.s, function(tau) vapply(copObs, iTau, numeric(1), tau = tau))
) # tTau is printed in moments.R

suppressWarnings(RNGversion("3.5.0")) ; set.seed(12)
u <- matrix(mku(1000), ncol= 2)# d = 2 required in the copula objects

##' == function(u) copula:::outside.01(u, strictly = FALSE) --- used in dCopula()
u.outside.01 <- function(u) apply(u, 1, function(x) any(x <= 0, 1 <= x))
u.out <- u.outside.01(u) # has NAs
u.ina <- apply(u, 1, anyNA)
u.iNA <- is.na(u.out)
table(u.ina, u.out, exclude={})
##        u.out
## u.ina   FALSE TRUE <NA>
##   FALSE   384   67    0
##   TRUE      0    7   42
u.OUT <- u.out & !u.ina  # no  NAs
## NB:  f(u) = 0 if some u_j is outside (0,1) - even when another u_k is NA

copsT <- lapply(setNames(,names(copObs)), # so have *named* list
                 function(cN)
                     lapply(unique(tTau[cN,]), function(th) setPar(copObs[[cN]], th)))
copsTh <- unlist(copsT, recursive=FALSE)
## --> a list of 72 parametrized copulas
okC <- lapply(copsTh, validObject, test=TRUE)
(inOk <- which(notOk <- !sapply(okC, isTRUE)))
## joeCopula2
##         43
if(length(inOk)) { # maybe no longer in future
    stopifnot(identical(inOk, c(joeCopula2 = 43L)))
    jp <- copsTh[[inOk]]@parameters
    print(jp - 1) # -3.41e-10
    stopifnot(all.equal(jp, 1))
    ## now fix it up:
    copsTh[[inOk]]@parameters <- 1
    validObject(copsTh[[inOk]])
}



## now using  ( u[], u.ina, u.OUT ):
op <- options(warn = 1) # immediate
## ==>  *all* Huessler-Reiss  produce warnings  such as  << FIXME eventually
##  Warning in log(u1p/u2p) : NaNs produced
## all other copulas give no warning !
errH <- if(getRversion() < "3.5.0") identity else function(e) {e$call <- NULL ;  e}
copChk <- lapply(copsTh, function(cop) {
    show(cop)
    fu <- dCopula(u, cop)
    Fu <- pCopula(u, cop)
    lfu <- dCopula(u, cop, log=TRUE)
    ## lfu. <- lfu[!u.ina]
    cat("-----------\n")
    tryCatch(
        stopifnot(fu[u.OUT] == 0,
		  all.equal(fu, exp(lfu), tolerance = 1e-15),
                  0 <= fu[!u.ina],
                  0 <= Fu[!u.ina], Fu[!u.ina] <= 1,
                  is.finite(lfu[!u.ina]) | lfu[!u.ina] == -Inf,
                  is.na(fu[u.iNA]), is.na(Fu[u.ina]), is.na(lfu[u.iNA]))
       , error = errH)
    })
options(op)

### FIXME:  Should see nothing below ( <==> length(ccc) == 0) ! ------------------------------

ccc <- unlist(copChk) # ==> drop all those that have NULL, i.e. *no* error above :
names(ccc) <- sub(".message$", "", names(ccc))
structure(as.list(ccc), class="simple.list")

stopifnot(length(ccc) <= 17)

## claytonCopula2     is.na(Fu[u.ina]) are not all TRUE  ((indepCopula case; progress already))
## frankCopula2       0 <= Fu[!u.ina] are not all TRUE   ((indepCopula case; fix?

## galambosCopula1    is.na(fu[u.iNA]) are not all TRUE  ((indepCopula case; progress already))
## galambosCopula2    0 <= Fu[!u.ina] are not all TRUE
## galambosCopula3    0 <= Fu[!u.ina] are not all TRUE
## galambosCopula4    0 <= Fu[!u.ina] are not all TRUE
## galambosCopula5    0 <= Fu[!u.ina] are not all TRUE

## huslerReissCopula1 0 <= Fu[!u.ina] are not all TRUE
## huslerReissCopula2 0 <= Fu[!u.ina] are not all TRUE
## huslerReissCopula3 0 <= Fu[!u.ina] are not all TRUE
## huslerReissCopula4 0 <= Fu[!u.ina] are not all TRUE
## huslerReissCopula5 0 <= Fu[!u.ina] are not all TRUE

## tawnCopula1        is.na(Fu[u.ina]) are not all TRUE
## tawnCopula2        is.na(Fu[u.ina]) are not all TRUE
## tawnCopula3        is.na(Fu[u.ina]) are not all TRUE
## tawnCopula4        is.na(Fu[u.ina]) are not all TRUE
## tawnCopula5        is.na(Fu[u.ina]) are not all TRUE

###-------------- "Unit" Tests for specific cases: ------------------------------

n.c <- as.integer(if(doExtras) 1001 else 101)

##== plackettCopula:  density log(dCopula()) vs "better"  dCopula(*, log=TRUE) ----- _UNFINISHED_ TODO

##' difference   log(c(..)) -  c(.., log=TRUE)
mkDfn <- function(u) {
    stopifnot(is.numeric(u), length(u) == 2 || (is.matrix(u) && ncol(u) == 2))
    Vectorize(function(alp) { C <- plackettCopula(alp); log(dCopula(u,C)) - dCopula(u,C,log=TRUE) })
}
Df <- mkDfn(u = c(.01,.99))
curve(Df, 0.99, 1.01, n=n.c) # alpha ~= 1 -- log=TRUE should be better -- but difference seem random!
curve(Df, 0,    2,    n=n.c) # to my surprise, alpha ~= 0 makes big difference
curve(Df, 1e-9, .01, log="x")
curve(abs(Df(x)), 1e-14, .1, log="xy")# and which one is more accurate ??
# another u[]
Df <- mkDfn(u = c(.1,.2))
curve(Df, 0.99, 1.01, n=n.c) # alpha ~= 1 -- log=TRUE: diff. larger for 1+eps than 1-eps
curve(Df, 0,    2,    n=n.c) # all very small

# and another
curve(mkDfn(u = c(.01,.02 ))(x), 0,  2,   n=n.c) # all very small
curve(mkDfn(u = c(.98,.999))(x), 0,  2,   n=n.c) # growing with alpha
curve(mkDfn(u = c(.98,.999))(x), .1, 1e4, n=n.c, log="x") # growing with alpha, still |.| < 6e-12
curve(mkDfn(u = c(.01,.001))(x), .01,1e8, n=n.c, log="x") # (growing) still small: |.| < 6e-15

##' 2-dim u[] in the low-density corner
mku <- function(n, one. = 0.999, eps = c(4^-(2:n2), b.n^-((n-n2):n))) {
    ## find basis for eps for very small eps
    stopifnot(n >= 2)
    n2 <- n %/% 2
    R.n <- - .Machine$double.min.exp / n
    b.n <- trunc(100*(2^R.n + 0.05))/100
    stopifnot(.5 > eps, eps > 0, 0.9 <= one., one. <= 1)
    cbind(one., eps, deparse.level = 0L)
}
if(doExtras) withAutoprint({
    alpha <- 1e6 ## <- large alpha
    (cop <- plackettCopula(alpha))
    summary(u <- mku(64))
    d1 <- dCopula(u, cop, log=TRUE)
    d2 <- log(dCopula(u, cop))
    summary(d1-d2)
    all.equal(d1,d2, tol=0)
})

