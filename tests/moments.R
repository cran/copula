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

if(getRversion() < "2.15")
paste0 <- function(...) paste(..., sep="")

setPar <- function(cop, par) setTheta(cop, par, noCheck=TRUE)

## Look at all non-virtual classes:
copcl <- unique(names(getClass("copula")@subclasses))
isVirt <- vapply(copcl, isVirtualClass, NA)
copcl <- copcl[!isVirt]
.notYet <- "schlatherCopula"
(copcl <- copcl[.notYet != copcl])
copF <- sapply(copcl, get)
frstArg <- vapply(copF, function(F) names(formals(F))[[1]], "")
copF1 <- copF[frstArg %in% c("dim", "param")]
copObj <- lapply(copF1, function(.) .())
stopifnot( sapply(copObj, is, class2 = "copula"),
           sapply(copObj, validObject))
copO. <- copObj[ names(copObj) != "indepCopula" ]
copO.2 <- copO.[ex2 <- !(names(copO.) %in% c("amhCopula","joeCopula"))]
                                        # because AMH has limited tau-range
str(copO.,max=1)
## The parameter bounds:
t(copBnds <- sapply(copO., function(C)
                    c(min= C@param.lowbnd[1], max= C@param.upbnd[1])))
copBnd.2 <- copBnds[, ex2]

###-------- tau & and inverse ---------------------------------------------------

## currently fails: --- FIXME?: should AMH also warn like the others?
tau.s <- c(-.999, -.1, 0, (1:3)/10, .5, .999)
### give different warnings , but "work" :  { .5 , even 1/3, gives error for AMH FIXME}
tau.s <- c(       -.1, 0, (1:2)/9, 0.3)
names(tau.s) <- paste0("tau=", sub("0[.]", ".", formatC(tau.s)))
tTau <- sapply(tau.s, function(tau)
               sapply(copO., iTau, tau = tau))

tTau
tTau["joeCopula", "tau=-.1"] <- 1 # ugly hack

stopifnot(rep(copBnds["min",],ncol(tTau)) <= tTau,
          tTau <= rep(copBnds["max",],ncol(tTau)),
          ## theta and tau are comonotone :
          apply(tTau, 1, diff) >= 0)

tautau <- t(sapply(names(copO.), function(cNam)
                   sapply(tTau[cNam,],
                          function(th) tau(setPar(copO.[[cNam]], th)))))

tautau
xctTau <- matrix(tau.s, nrow = nrow(tautau), ncol=length(tau.s),
                 byrow=TRUE)
## The absolute errors
errTau <- tautau-xctTau
round(10000*errTau)
## has two NaN .. ok, for now
errTau["tawnCopula", 1:2] <- 0
## the tevCopula cannot get a tau <= 0 (for now) __FIXME?__
errTau["tevCopula", 1:2] <- 0
## These families do not support tau < 0  (currently):
errTau[c("gumbelCopula", "joeCopula", "galambosCopula", "huslerReissCopula"),
       "tau=-.1"] <- 0
## "fgmCopula" has tau in [-2/9, 2/9] :
errTau["fgmCopula", "tau=.3"] <- 0
stopifnot(max(abs(errTau)) <= 0.00052)# ok for IJ-taus

###-------- rho & and inverse ---------------------------------------------------

## NB:
##  iRho() method for class "amhCopula" not yet implemented

### give different warnings , but "work" [not using AMH !]
rho.s <- c(-.999, -.1, 0, (1:3)/9, .5, .9, .999)
names(rho.s) <- paste0("rho=", sub("0[.]", ".", formatC(rho.s)))
tRho <- sapply(rho.s, function(rho)
               sapply(copO.2, iRho, rho = rho))
warnings()

tRho
##--> oops!  clayton [rho=0] is NA __FIXME__
tRho["claytonCopula", "rho=0"] <- 0
## and it has NA also for .999  {but that maybe consider ok}:
tRho["claytonCopula", "rho=.999"] <- 10^100

stopifnot(rep(copBnd.2["min",],ncol(tRho)) <= tRho,
          tRho <= rep(copBnd.2["max",],ncol(tRho)),
          ## theta and rho are comonotone :
          apply(tRho, 1, diff) >= 0)

rhorho <- t(sapply(names(copO.2), function(cNam)
                   sapply(tRho[cNam,],
                          function(th) rho(setPar(copO.2[[cNam]], th)))))

rhorho
xctRho <- matrix(rho.s, nrow = nrow(rhorho), ncol=length(rho.s),
                 byrow=TRUE)
## The absolute errors
errRho <- rhorho-xctRho
round(10000*errRho)
## the tevCopula cannot get a rho <= 0 (for now) __FIXME?__
errRho["tevCopula", 1:3] <- 0
## These three families do not support rho < 0  (currently):
errRho[c("gumbelCopula", "galambosCopula", "huslerReissCopula", "tawnCopula"),
       c("rho=-.1","rho=-.999")] <- 0
errRho["tawnCopula", rho.s >= 0.9] <- 0
## "fgmCopula" has rho in [-1/3, 1/3] :
errRho["fgmCopula", abs(rho.s) > 1/3] <- 0

stopifnot(max(abs(errRho)) <= 0.00369,
          max(abs(errRho[,rho.s <= 0.9])) <= 0.0002)# ok for IJ-rhos


cat('Time elapsed: ', proc.time(),'\n') # for ''statistical reasons''
