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


### Quasi-random numbers for copula models #####################################

## Note: see ``A primer on quasi-random numbers for copula models'' for more
##       details.


### 0) Setup ###################################################################

require(lattice)
require(copula)
require(VineCopula)
require(randtoolbox) # for quasi-random number generators
require(VGAM)

doPDF <- !dev.interactive(orNone=TRUE)


### 1) Functions ###############################################################

##' @title Inverse of the bivariate conditional Marshall--Olkin copula
##' @param u (n,2) matrix of U[0,1] random numbers to be transformed to
##'        (u[,1], C^-(u[,2]|u[,1]))
##' @param alpha bivariate parameter vector
##' @return (u[,1], C^-(u[,2]|u[,1])) for C being a MO copula
##' @author Marius Hofert
C.inv.MO <- function(u, alpha)
{
    stopifnot(is.matrix(u), 0 <= alpha, alpha <= 1)
    up <- u[,1]^(alpha[1]*(1/alpha[2]-1))
    low <- (1-alpha[1])*up
    i1 <- u[,2] <= low
    i3 <- u[,2] >  up
    u2 <- u[,1]^(alpha[1]/alpha[2])
    u2[i1] <- u[i1,1]^alpha[1] * u[i1,2] / (1-alpha[1])
    u2[i3] <- u[i3,2]^(1/(1-alpha[2]))
    cbind(u[,1], u2)
}

##' @title Sampling bivariate Clayton copulas via MO based on transformed U[0,1]^3
##' @param u (n,3) matrix of U[0,1] random numbers to be transformed
##' @param theta Clayton parameter
##' @return (n,2) matrix of transformed random numbers following a bivariate
##'         Clayton copula
##' @author Marius Hofert
MOtrafoC <- function(u, theta)
    copClayton@psi(-log(u[,1:2]) / qgamma(u[,3], 1/th), theta=th)

##' @title Constructing points on a square for plotting purposes
##' @param ll lower left endpoint of the square
##' @param ur upper right endpoint of the square
##' @param theta Clayton parameter
##' @param np number of points on the square in each of the four directions
##' @param col colors for the four sides of the squares
##' @return a list including
##'         u:     (4x np) points on the square
##'         label: a label indicating the side of the square
##'         col:   colors for the points on the square
##'         u.CDM: transformed square (via CDM)
##'         u.VWS: transformed square (via VWS)
##'         u.MO:  transformed square (via MO)
##' @author Marius Hofert
square <- function(ll, ur, theta, np=16,
                   col=c("darkorange2", "firebrick", "royalblue3", "black"))
{
    stopifnot(length(ll) == 2, length(ur) == 2, length(col) == 4,
              ll[1] < ur[1], ll[2] < ur[2], np > 1)
    ## build points
    SW.to.NW <- cbind(u1=rep(ll[1], np), u2=seq(ll[2], ur[2], length.out=np))
    NW.to.NE <- cbind(u1=seq(ll[1], ur[1], length.out=np), u2=rep(ur[2], np))
    NE.to.SE <- cbind(u1=rep(ur[1], np), u2=seq(ur[2], ll[2], length.out=np))
    SE.to.SW <- cbind(u1=seq(ur[1], ll[1], length.out=np), u2=rep(ll[2], np))
    u <- rbind(SW.to.NW, NW.to.NE, NE.to.SE, SE.to.SW)
    ## build attributes
    label <- rep(c("SW.to.NW", "NW.to.NE", "NE.to.SE", "SE.to.SW"), each=np)
    col <- rep(col, each=np)
    ## build transformed points
    cop <- onacopulaL("Clayton", nacList=list(theta, 1:2))
    u.CDM <- rtrafo(u, cop=cop, inverse=TRUE)
    u.WVS <- htrafo(u, cop=cop, inverse=TRUE)
    u.qrn <- halton(4*np) # additional sequence (deterministic!) for MO trafo
    u.MO  <- MOtrafoC(cbind(u, u.qrn), theta=theta)
    ## return
    list(u=u, label=label, col=col, u.CDM=u.CDM, u.WVS=u.WVS, u.MO=u.MO)
}


### 2) PRNG vs QRNG ############################################################

n <- 1000 # sample size


### Independence copula ########################################################

## PRNG
set.seed(271)
U <- matrix(runif(n*2), ncol=2)
if(doPDF) pdf(file=(file <- "fig_prng.pdf"), width=6, height=6)
par(pty="s")
plot(U, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)

## QRNG
set.seed(271)
U. <- halton(n, dim=2)
if(doPDF) pdf(file=(file <- "fig_qrng.pdf"), width=6, height=6)
par(pty="s")
plot(U., xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)


### t_3 copula #################################################################

family <- "t"
nu <- 3 # degrees of freedom
tau <- 0.5
th <- iTau(ellipCopula(family, df=nu), tau)
cop <- ellipCopula(family, param=th, df=nu)

## PRNG
U.t <- rtrafo(U, cop=cop, inverse=TRUE)
file <- paste0("fig_prng_t", nu, "_", tau, ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.t, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)

## QRNG
U.t. <- rtrafo(U., cop=cop, inverse=TRUE)
file <- paste0("fig_qrng_t", nu, "_", tau, ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.t., xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)

## QRNG (pairs plot for 3d t_3 sample)
set.seed(271)
U.3d. <- halton(n, dim=3)
cop3d <- ellipCopula(family, param=th, dim=3, df=nu)
U.t.3d. <- rtrafo(U.3d., cop=cop3d, inverse=TRUE)
file <- paste0("fig_qrng_t", nu, "_", tau, "_pairs.pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
pairs(U.t.3d., gap=0,
      labels=as.expression(sapply(1:3, function(j) bquote(italic(U[.(j)])))))
if(doPDF) dev.off.pdf(file=file)
## => projections (here: to pairs) can be deceiving/non-optimal
##    (but 'quasi-randomness' not easily visible from a 3d cloud plot either)


### Clayton copula #############################################################

family <- "Clayton"
tau <- 0.5
th <- iTau(getAcop(family), tau)
cop <- onacopulaL(family, nacList=list(th, 1:2))

## PRNG
U.C <- rtrafo(U, cop=cop, inverse=TRUE)
file <- paste0("fig_prng_C_", tau, ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.C, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)

## QRNG
U.C. <- rtrafo(U., cop=cop, inverse=TRUE)
file <- paste0("fig_qrng_C_", tau, ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.C., xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)


### Marshall--Olkin copula #####################################################

alpha <- c(0.25, 0.75)
tau <- (alpha[1]*alpha[2]) / (alpha[1]+alpha[2]-alpha[1]*alpha[2])

## PRNG
U.MO <- C.inv.MO(U, alpha=alpha)
file <- paste0("fig_prng_MO_", paste0(alpha, collapse="_"), ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.MO, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)

## QRNG
U.MO. <- C.inv.MO(U., alpha=alpha)
file <- paste0("fig_qrng_MO_", paste0(alpha, collapse="_"), ".pdf")
if(doPDF) pdf(file=file, width=6, height=6)
par(pty="s")
plot(U.MO., xlab=expression(italic(U[1])), ylab=expression(italic(U[2])))
if(doPDF) dev.off.pdf(file=file)


### R-Vine copula ##############################################################

## Note: RVineSim() [in RVineSim.R] -> SimulateRVine() [in rvine.c] ->
##       Hinv1() [in hfunc.c] -> Hinv()
##       => the following families require *numerical* root finding
##          (via HNumInv(); may destroy low discrepancy)
##          for computing the inverses of the h functions:
##          0: Pi (fine)
##          1: Ga (uses pnorm/qnorm -- should be fine, too)
##          2: t  (should be fine, too)
##          3: C  (fine)
##          4: G  (numerical inversion!)
##          5: F  (should be fine)
##          6: J  (numerical inversion!)


### 3d example #################################################################

## R-vine tree structure matrix
M <- matrix(c(3, 1, 2,
              0, 2, 1,
              0, 0, 1), ncol=3)

## R-vine pair-copula family matrix (0 = Pi)
family <- matrix(c(0, 3, 3, # C, C
                   0, 0, 3, #    C
                   0, 0, 0), ncol=3)

## define R-vine pair-copula parameter matrix
param <- matrix(c(0, 1, 1,
                  0, 0, 1,
                  0, 0, 0), ncol=3)

## define second R-vine pair-copula parameter matrix
param. <- matrix(0, nrow=3, ncol=3)

## define RVineMatrix object
RVM <- RVineMatrix(Matrix=M, family=family, par=param, par2=param.)

## PRNG
set.seed(271)
U <- RVineSim(n, RVM)
if(doPDF) pdf(file=(file <- "fig_R-vine_prng_d=3.pdf"), width=6, height=6)
pairs(U, labels=as.expression( sapply(1:3, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)

## QRNG (Halton)
set.seed(271)
U. <- halton(n, d=3)
if(doPDF) pdf(file=(file <- "fig_qrng_d=3_halton.pdf"), width=6, height=6)
pairs(U., labels=as.expression( sapply(1:3, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)
## transform to copula data
U. <- RVineSim(n, RVM, U=U.)
if(doPDF) pdf(file=(file <- "fig_R-vine_qrng_d=3_halton.pdf"), width=6, height=6)
pairs(U., labels=as.expression( sapply(1:3, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)

## QRNG (Sobol)
set.seed(271)
U. <- sobol(n, d=3)
if(doPDF) pdf(file=(file <- "fig_qrng_d=3_sobol.pdf"), width=6, height=6)
pairs(U., labels=as.expression( sapply(1:3, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)
## transform to copula data
U. <- RVineSim(n, RVM, U=U.)
if(doPDF) pdf(file=(file <- "fig_R-vine_qrng_d=3_sobol.pdf"), width=6, height=6)
pairs(U., labels=as.expression( sapply(1:3, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)
## Similarly to the 3d t copula case, by the 'projection to pairs argument'
## not all pairs look 'quasi-random' (note: for some family choices, there is
## numerical root-finding involved, but not in this case here)


### 5d example #################################################################

## R-vine tree structure matrix
M <- matrix(c(5, 2, 3, 1, 4,
              0, 2, 3, 4, 1,
              0, 0, 3, 4, 1,
              0, 0, 0, 4, 1,
              0, 0, 0, 0, 1), ncol=5)

## R-vine pair-copula family matrix (0 = Pi)
family <- matrix(c(0, 1, 3, 4, 4, # Ga, C, G, G
                   0, 0, 3, 4, 1, #     C, G, Ga
                   0, 0, 0, 4, 1, #        G, Ga
                   0, 0, 0, 0, 3, #           C
                   0, 0, 0, 0, 0), ncol=5)

## define R-vine pair-copula parameter matrix
param <- matrix(c(0, 0.2, 0.9, 1.5, 3.9,
                  0,   0, 1.1, 1.6, 0.9,
                  0,   0,   0, 1.9, 0.5,
                  0,   0,   0,   0, 4.8,
                  0,   0,   0,   0,   0), ncol=5)

## define second R-vine pair-copula parameter matrix
param. <- matrix(0, nrow=5, ncol=5)

## define RVineMatrix object
RVM <- RVineMatrix(Matrix=M, family=family, par=param, par2=param.)

## PRNG
set.seed(271)
U <- RVineSim(n, RVM)
if(doPDF) pdf(file=(file <- "fig_R-vine_prng_d=5.pdf"), width=6, height=6)
pairs(U, labels=as.expression( sapply(1:5, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)

## QRNG
set.seed(271)
U. <- RVineSim(n, RVM, U=halton(n, d=5))
if(doPDF) pdf(file=(file <- "fig_R-vine_qrng_d=5.pdf"), width=6, height=6)
pairs(U., labels=as.expression( sapply(1:5, function(j) bquote(italic(U[.(j)]))) ),
      gap=0, cex=0.3)
if(doPDF) dev.off.pdf(file=file)


### 3) Why non-one-to-one transformations (may) fail ###########################

## parameters declaration
n <- 1000
family <- "Clayton"
tau <- 0.5
th <- iTau(getAcop(family), tau)
cop <- onacopulaL(family, nacList=list(th, 1:2))


### 3.1) Colorized scatter plot ################################################

## sample quasi-random numbers for CDM, WVS
set.seed(271)
U.. <- halton(n, 3) # 3 due to MO
U.  <- U..[,1:2] # for CDM, WVS to be as 'close/comparable' to MO as possible

## assign corresponding colors
col <- rep("black", n)
col[U.[,1] <= 0.4 & U.[,2] <= 0.4] <- "firebrick"
col[U.[,1] >= 0.8 & U.[,2] >= 0.8] <- "royalblue3"

## transform into dependent samples
U_CDM <- rtrafo(U., cop=cop, inverse=TRUE)
U_WVS <- htrafo(U., cop=cop, inverse=TRUE)
U_MO  <- MOtrafoC(U.., theta=th)

## Colorized scatter plot (quasi-random numbers)
if(doPDF) pdf(file=(file <- "fig_qrng_col.pdf"), width=6, height=6)
par(pty="s")
plot(U., xlab=expression(italic(U[1])), ylab=expression(italic(U[2])), col=col)
if(doPDF) dev.off.pdf(file=file)

## Colorized scatter plots (CDM)
if(doPDF) pdf(file=(file <- "fig_qrng_col_CDM.pdf"), width=6, height=6)
par(pty="s")
plot(U_CDM, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])), col=col)
if(doPDF) dev.off.pdf(file=file)

## Colorized scatter plots (WVS)
if(doPDF) pdf(file=(file <- "fig_qrng_col_WVS.pdf"), width=6, height=6)
par(pty="s")
plot(U_WVS, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])), col=col)
if(doPDF) dev.off.pdf(file=file)

## Colorized scatter plots (MO)
if(doPDF) pdf(file=(file <- "fig_qrng_col_MO.pdf"), width=6, height=6)
par(pty="s")
plot(U_MO, xlab=expression(italic(U[1])), ylab=expression(italic(U[2])), col=col)
if(doPDF) dev.off.pdf(file=file)


### 3.2) Mapping nested squares ################################################

## construct squares
sq.out <- square(ll=c(0.2,0.2), ur=c(0.8,0.8), theta=th)
sq.mid <- square(ll=c(0.3,0.3), ur=c(0.7,0.7), theta=th)
sq.in  <- square(ll=c(0.4,0.4), ur=c(0.6,0.6), theta=th)

## Squares
if(doPDF) pdf(file=(file <- "fig_squares.pdf"), width=6, height=6)
par(pty="s")
plot(NULL, type="l", xlim=0:1, ylim=0:1,
     xlab=expression(italic(u[1])), ylab=expression(italic(u[2])))
for(i in 1:(nrow(sq.out[["u"]])-1)) lines(sq.out[["u"]][i:(i+1),], col=sq.out[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.mid[["u"]])-1)) lines(sq.mid[["u"]][i:(i+1),], col=sq.mid[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.in[["u"]])-1)) lines(sq.in[["u"]][i:(i+1),], col=sq.in[["col"]][i], lwd=2)
if(doPDF) dev.off.pdf(file=file)

## CDM-transformed squares
if(doPDF) pdf(file=(file <- "fig_squares_CDM.pdf"), width=6, height=6)
par(pty="s")
plot(NULL, type="l", xlim=0:1, ylim=0:1,
     xlab=expression(italic(u[1])), ylab=expression(italic(u[2])))
for(i in 1:(nrow(sq.out[["u.CDM"]])-1)) lines(sq.out[["u.CDM"]][i:(i+1),], col=sq.out[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.mid[["u.CDM"]])-1)) lines(sq.mid[["u.CDM"]][i:(i+1),], col=sq.mid[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.in[["u.CDM"]])-1)) lines(sq.in[["u.CDM"]][i:(i+1),], col=sq.in[["col"]][i], lwd=2)
if(doPDF) dev.off.pdf(file=file)

## WVS-transformed squares
if(doPDF) pdf(file=(file <- "fig_squares_WVS.pdf"), width=6, height=6)
par(pty="s")
plot(NULL, type="l", xlim=0:1, ylim=0:1,
     xlab=expression(italic(u[1])), ylab=expression(italic(u[2])))
for(i in 1:(nrow(sq.out[["u.WVS"]])-1)) lines(sq.out[["u.WVS"]][i:(i+1),], col=sq.out[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.mid[["u.WVS"]])-1)) lines(sq.mid[["u.WVS"]][i:(i+1),], col=sq.mid[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.in[["u.WVS"]])-1)) lines(sq.in[["u.WVS"]][i:(i+1),], col=sq.in[["col"]][i], lwd=2)
if(doPDF) dev.off.pdf(file=file)

## MO-transformed squares
if(doPDF) pdf(file=(file <- "fig_squares_MO.pdf"), width=6, height=6)
par(pty="s")
plot(NULL, type="l", xlim=0:1, ylim=0:1,
     xlab=expression(italic(u[1])), ylab=expression(italic(u[2])))
for(i in 1:(nrow(sq.out[["u.MO"]])-1)) lines(sq.out[["u.MO"]][i:(i+1),], col=sq.out[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.mid[["u.MO"]])-1)) lines(sq.mid[["u.MO"]][i:(i+1),], col=sq.mid[["col"]][i], lwd=2)
for(i in 1:(nrow(sq.in[["u.MO"]])-1)) lines(sq.in[["u.MO"]][i:(i+1),], col=sq.in[["col"]][i], lwd=2)
if(doPDF) dev.off.pdf(file=file)

