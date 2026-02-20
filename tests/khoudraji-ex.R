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

### --- Khoudraji Copula -----------------------

library(copula)

source(system.file("Rsource", "utils.R", package="copula", mustWork=TRUE))
## tryCatch.W.E() etc

(doExtras <- copula:::doExtras())

isExplicit <- copula:::isExplicit
## if(!dev.interactive(orNone=TRUE)) pdf("khoudraji-ex.pdf")

################################################################################
## An explicit Khoudraji copula constructed from claytonCopula and indepCoula
## with one fixed shape parameter
kcf <- khoudrajiCopula(copula2 = claytonCopula(6),
                       shapes = fixParam(c(0.4, 0.95), c(FALSE, TRUE)))
kcf
stopifnot(
    all.equal(getTheta(kcf, freeOnly = FALSE, named = TRUE) -> th.,
              c(c2.alpha = 6, shape1 = 0.4, shape2 = 0.95)),
    identical(getTheta(kcf, named = TRUE), th.[1:2]))

## kcf@exprdist$cdf
stopifnot(isExplicit(kcf))

set.seed(123)
u <- rCopula(20, kcf)

## cdf/pdf from explicit expression versus from algorithmic implementation
summary(rD  <- pCopula(u, kcf) - copula:::pKhoudrajiCopula(u, kcf))
summary(rD. <- dCopula(u, kcf) - copula:::dKhoudrajiBivCopula(u, kcf))
stopifnot(max(abs(rD)) < 1e-12, max(abs(rD.)) < 1e-12)
showProc.time()

################################################################################
## A nested Khoudraji copula from kcf and a gumbelCopula, still explicit
k_kcf_g <- khoudrajiCopula(kcf, gumbelCopula(2),
                           shapes = fixParam(c(.2, .9), c(TRUE, TRUE)))
k_kcf_g
stopifnot(isExplicit(k_kcf_g))
getTheta(k_kcf_g, freeOnly = FALSE, named = TRUE)
getTheta(k_kcf_g, named = TRUE)

## k_kcf_g@exprdist$cdf
set.seed(456)
U <- rCopula(100000, k_kcf_g)

## check cdf/pdf with C.n and kde2d
require(MASS)
u <- as.matrix(expand.grid((1:9)/10, (1:9)/10))
## cdf versus C.n
summary(dp <- pCopula(u, k_kcf_g) - C.n(u, U)) ## |.| <= 9.8e-4)
stopifnot(max(abs(dp)) < 3e-3)
## pdf versus kde2d
kde <- kde2d(U[,1], U[,2], n = 9, lims = c(0.1, 0.9, 0.1, 0.9))
summary(rD <- dCopula(u, k_kcf_g) / c(kde$z) - 1) ## |relative difference| <= 0.123..
stopifnot(max(abs(rD)) < 0.13)
showProc.time()

################################################################################
## one more nesting: kcf and k_kcf_g
monster <- khoudrajiCopula(kcf, k_kcf_g,
                           shapes = fixParam(c(.9, .1), c(TRUE, FALSE)))
monster; validObject(monster)
stopifnot(isExplicit(monster))

(th.mA <- getTheta(monster, freeOnly = FALSE, named = TRUE))
(th.m <- getTheta(monster, named = TRUE))
stopifnot(identical(th.m, th.mA[isFree(monster)]))

set.seed(488)
U <- rCopula(100000, monster)
## cdf versus C.n
summary(dp <- pCopula(u, monster) - C.n(u, U)) ## -0.0083 ... 0.00131..
## pdf versus kde2d
kde <- kde2d(U[,1], U[,2], n = 9, lims = c(0.1, 0.9, 0.1, 0.9))
summary(rdC <- dCopula(u, monster) / c(kde$z) - 1) ## relative difference < 0.09
stopifnot(abs(dp) < 0.01, abs(rdC) < 0.09)
showProc.time()

################################################################################
## khoudrajiExplicitCopula with dim 3
kcd3 <- khoudrajiCopula(copula1 = indepCopula(dim=3),
                        copula2 = claytonCopula(6, dim=3),
                        shapes = c(0.4, 0.95, 0.95))
kcd3

set.seed(1712)
U <- rCopula(100000, kcd3)
u <- matrix(runif(15), 5, 3)
## cdf versus C.n
summary(dP <- pCopula(u, kcd3) - C.n(u, U)) # < 0.00007
## don't know how to check pdf
(f.v <- dCopula(u, kcd3))
stopifnot(expr = {
    isExplicit(kcd3)
    inherits(kcd3, "khoudrajiExplicitCopula")
    all.equal(f.v,   tolerance = 1e-7,
              c(1.1116773, 1.6661734, 0.1719080, 0.4981574, 9.8964415))
    abs(dP) < 8e-4
})
showProc.time()

kcg <- khoudrajiCopula(claytonCopula(2, dim=3), gumbelCopula(3, dim=3), c(0.2,  0.2, .8))
kcg
u <- rCopula(10, kcg)
dCopula(u, kcg)


kkcgg <- khoudrajiCopula(kcg, gumbelCopula(2, dim = 3), c(.3, .3, .9))
kkcgg
u <- rCopula(10, kkcgg)
dCopula(u, kkcgg)


kcgcd3 <- khoudrajiCopula(kcd3, kkcgg, c(.4, .4, .5))
kcgcd3
u <- rCopula(10, kcgcd3)
dCopula(u, kcgcd3)

showProc.time()
if(doExtras) {
 print(compD <- comparederiv(kcgcd3, u), digits = 4)
 ##         dCdu     dCdtheta      dlogcdu  dlogcdtheta
 ## 4.387280e-03 9.791440e-04 3.009563e-08 7.572446e-09
 stopifnot(compD < c(0.01, 3e-3, 1e-7, 8e-8))
 showProc.time()
}
#############################################################################
## fitting checking

do1 <- function(n, cop) {
    U <- pobs(rCopula(n, cop))
    fit <- try(fitCopula(cop, U, method = "mpl"))
    p <- nParam(cop, freeOnly=TRUE)
    if (inherits(fit, "try-error")) rep(NA_real_, 2 * p)
    else c(coef(fit), sqrt(diag(vcov(fit))))
}

sumsim <- function(sim) {
    p <- nrow(sim) / 2
    mm <- apply(sim, 1, mean, na.rm = TRUE)
    ss <- apply(sim, 1, sd, na.rm = TRUE)
    data.frame(estimate = mm[1:p], ese = ss[1:p], ase = mm[p + 1:p])
}

if (doExtras) {
    set.seed(123)
    nrep <- 50
    n <- 1000

    sim <- replicate(nrep, do1(n, kcg))
    sumsim(sim)
    cbind(true = getTheta(kcg), sumsim(sim))
    ## the average se is greater than empirical for the shape parameters; need more investigation
}
